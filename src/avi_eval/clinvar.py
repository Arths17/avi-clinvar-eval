"""Read a ClinVar GRCh38 VCF and apply the pre-registered filters.

Plain-English summary of what happens here:
  1. Stream the VCF line by line (it is large) and count every line.
  2. Keep single-nucleotide variants (SNVs) on chromosomes 1-22 and X.
  3. Keep only variants whose ClinVar classification is cleanly pathogenic-type
     or cleanly benign-type (no uncertain / conflicting / risk-factor records).
  4. Convert ClinVar's review-status text to a 0-4 star rating.
  5. Assign each variant a consequence category ("stratum") from its MC field.
Every step's remaining count is written to an attrition log.
"""
from __future__ import annotations

import gzip
import logging
import re
from collections import Counter
from dataclasses import dataclass, field
from pathlib import Path
from typing import Iterator

import pandas as pd

from . import config

log = logging.getLogger(__name__)

REQUIRED_INFO_IDS = ("CLNSIG", "CLNREVSTAT", "MC", "GENEINFO", "CLNVC")
_BASES = frozenset("ACGT")
_INFO_ID_RE = re.compile(r"^##INFO=<ID=([^,>]+)")


# --------------------------------------------------------------------------- #
# Attrition log
# --------------------------------------------------------------------------- #
@dataclass
class Attrition:
    """Ordered record of how many variants remain after each filter step."""

    rows: list[dict] = field(default_factory=list)

    def add(self, step: str, n: int, note: str = "") -> None:
        prev = self.rows[-1]["n_remaining"] if self.rows else None
        self.rows.append(
            {"step": step, "n_remaining": int(n), "n_dropped": None if prev is None else int(prev - n), "note": note}
        )
        log.info("ATTRITION %-45s remaining=%d %s", step, n, note)

    def to_frame(self) -> pd.DataFrame:
        return pd.DataFrame(self.rows)

    def save(self, path: Path) -> None:
        path.parent.mkdir(parents=True, exist_ok=True)
        self.to_frame().to_csv(path, index=False)


# --------------------------------------------------------------------------- #
# Small pure helpers (unit-tested)
# --------------------------------------------------------------------------- #
def normalise_chrom(chrom: str) -> str:
    """'chr1' -> '1', '1' -> '1'."""
    return chrom[3:] if chrom.lower().startswith("chr") else chrom


def is_snv(ref: str, alt: str) -> bool:
    """True for a single reference base and a single, different alternate base (A/C/G/T)."""
    return len(ref) == 1 and len(alt) == 1 and ref in _BASES and alt in _BASES and ref != alt


def label_from_clnsig(clnsig: str | None) -> int | None:
    """Return 1 (pathogenic-type), 0 (benign-type) or None (excluded).

    The CLNSIG string may hold several terms separated by '|' or ','. Every term
    must belong to the same polarity set; anything else (uncertain, conflicting,
    risk factor, drug response, low penetrance, mixed polarity) is excluded.
    """
    if not clnsig:
        return None
    terms = [t.strip().lower() for t in re.split(r"[|,]", clnsig) if t.strip()]
    if not terms:
        return None
    if all(t in config.POSITIVE_TERMS for t in terms):
        return 1
    if all(t in config.NEGATIVE_TERMS for t in terms):
        return 0
    return None


def stars_from_revstat(revstat: str | None) -> int:
    """ClinVar review status text -> 0-4 stars; -1 if the text is not recognised."""
    if revstat is None:
        return -1
    key = revstat.strip().lower().replace(",", "")
    return config.REVIEW_STARS.get(key, -1)


def parse_mc(mc: str | None) -> list[str]:
    """'SO:0001583|missense_variant,SO:0001627|intron_variant' -> ['missense_variant', 'intron_variant']."""
    if not mc:
        return []
    out = []
    for item in mc.split(","):
        _, _, term = item.partition("|")
        term = (term or item).strip().lower()
        if term:
            out.append(term)
    return out


def assign_stratum(terms: list[str]) -> str | None:
    """First stratum (in priority order) whose term set intersects `terms`; None if unmapped."""
    tset = set(terms)
    for stratum, stratum_terms in config.STRATUM_TERMS.items():
        if tset & stratum_terms:
            return stratum
    return None


def first_gene(geneinfo: str | None) -> str | None:
    """'BRCA1:672|BRCA1-AS:1234' -> 'BRCA1' (design.md deviation 3)."""
    if not geneinfo:
        return None
    first = geneinfo.split("|")[0]
    symbol = first.split(":")[0].strip()
    return symbol or None


def check_header(header_lines: list[str]) -> None:
    """Fail loudly if the VCF does not declare the INFO fields this study depends on."""
    declared = {m.group(1) for line in header_lines if (m := _INFO_ID_RE.match(line))}
    missing = [i for i in REQUIRED_INFO_IDS if i not in declared]
    if missing:
        raise ValueError(
            f"ClinVar VCF header is missing INFO fields {missing}. "
            f"Declared fields: {sorted(declared)}. Stop and check the file/format before continuing."
        )


# --------------------------------------------------------------------------- #
# VCF streaming
# --------------------------------------------------------------------------- #
def _open_text(path: Path):
    return gzip.open(path, "rt", encoding="utf-8") if str(path).endswith(".gz") else open(path, encoding="utf-8")


def _pick_info(info: str, keys: tuple[str, ...]) -> dict[str, str]:
    out: dict[str, str] = {}
    for item in info.split(";"):
        k, sep, v = item.partition("=")
        if sep and k in keys:
            out[k] = v
    return out


def load_clinvar(vcf_path: Path, attrition: Attrition) -> tuple[pd.DataFrame, Counter]:
    """Stream the VCF, apply filters up to (and including) the star rating, and
    return a DataFrame with one row per labelled SNV with stars >= 1.

    Also returns a Counter of unrecognised review-status strings (so nothing is
    silently dropped).
    """
    header: list[str] = []
    n_raw = n_snv = n_chrom = n_label = 0
    unknown_revstat: Counter = Counter()
    kept: list[dict] = []
    keys = ("CLNSIG", "CLNREVSTAT", "MC", "GENEINFO")

    with _open_text(vcf_path) as fh:
        for line in fh:
            if line.startswith("##"):
                header.append(line.rstrip("\n"))
                continue
            if line.startswith("#"):
                check_header(header)
                continue
            n_raw += 1
            f = line.rstrip("\n").split("\t", 8)
            if len(f) < 8:
                continue
            chrom, pos, vid, ref, alt, info = normalise_chrom(f[0]), f[1], f[2], f[3], f[4], f[7]
            if not is_snv(ref, alt):
                continue
            n_snv += 1
            if chrom not in config.CHROMS:
                continue
            n_chrom += 1
            d = _pick_info(info, keys)
            label = label_from_clnsig(d.get("CLNSIG"))
            if label is None:
                continue
            n_label += 1
            stars = stars_from_revstat(d.get("CLNREVSTAT"))
            if stars < 0:
                unknown_revstat[d.get("CLNREVSTAT", "<missing>")] += 1
                continue
            if stars < 1:
                continue
            kept.append(
                {
                    "variation_id": int(vid),
                    "chrom": chrom,
                    "pos": int(pos),
                    "ref": ref,
                    "alt": alt,
                    "label": label,
                    "stars": stars,
                    "mc": d.get("MC", ""),
                    "gene": first_gene(d.get("GENEINFO")),
                }
            )
    if not header:
        raise ValueError(f"{vcf_path} has no VCF header lines")

    df = pd.DataFrame(kept)
    attrition.add("raw VCF data rows", n_raw)
    attrition.add("single-nucleotide variants", n_snv)
    attrition.add("on chr1-22 or X", n_chrom)
    attrition.add("clean pathogenic-type or benign-type label", n_label)
    attrition.add("recognised review status and >=1 star", len(df), note=f"unrecognised review status: {sum(unknown_revstat.values())}")
    return df, unknown_revstat


def add_strata(df: pd.DataFrame) -> pd.DataFrame:
    """Add `stratum` (None where unmapped) and `heldout_chrom` columns."""
    out = df.copy()
    out["stratum"] = [assign_stratum(parse_mc(m)) for m in out["mc"]]
    out["heldout_chrom"] = out["chrom"].isin(config.HELDOUT_CHROMS)
    return out


def mc_term_counts(df: pd.DataFrame) -> pd.DataFrame:
    """Counts of every MC term by class - used at the stratum-mapping gate (counts only, no scores)."""
    rows: Counter = Counter()
    for mc, label in zip(df["mc"], df["label"]):
        for term in set(parse_mc(mc)):
            rows[(term, int(label))] += 1
    out = pd.DataFrame([{"mc_term": t, "label": lab, "n": n} for (t, lab), n in rows.items()])
    if out.empty:
        return pd.DataFrame(columns=["mc_term", "n_pathogenic", "n_benign", "mapped_to"])
    wide = out.pivot_table(index="mc_term", columns="label", values="n", fill_value=0).reset_index()
    wide = wide.rename(columns={1: "n_pathogenic", 0: "n_benign"})
    for c in ("n_pathogenic", "n_benign"):
        if c not in wide:
            wide[c] = 0
    wide["mapped_to"] = [assign_stratum([t]) or "UNMAPPED" for t in wide["mc_term"]]
    return wide.sort_values("n_pathogenic", ascending=False)[["mc_term", "n_pathogenic", "n_benign", "mapped_to"]]


def iter_download_chunks(response, chunk: int = 1 << 20) -> Iterator[bytes]:
    while True:
        b = response.read(chunk)
        if not b:
            return
        yield b
