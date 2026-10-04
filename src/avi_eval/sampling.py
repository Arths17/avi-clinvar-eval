"""Reproducible stratified sampling (design.md section 4).

For every (stratum, class) cell:
  1. shuffle the eligible variants with a seeded generator,
  2. keep at most GENE_CAP variants per gene (variants with no gene are exempt),
  3. keep the first N of what remains (so the final pick is a uniform random subset).
"""
from __future__ import annotations

import numpy as np
import pandas as pd

from . import config
from .clinvar import Attrition

SAMPLE_COLUMNS = [
    "variation_id", "chrom", "pos", "ref", "alt", "label", "stars",
    "stratum", "gene", "heldout_chrom", "cohort",
]


def sample_cell(cell: pd.DataFrame, n: int, gene_cap: int, rng: np.random.Generator) -> pd.DataFrame:
    """Sample up to `n` rows from one (stratum, class) cell."""
    if cell.empty:
        return cell
    cell = cell.sort_values("variation_id").reset_index(drop=True)
    shuffled = cell.iloc[rng.permutation(len(cell))].reset_index(drop=True)
    has_gene = shuffled["gene"].notna()
    capped_with_gene = shuffled[has_gene].groupby("gene", sort=False).head(gene_cap)
    capped = pd.concat([capped_with_gene, shuffled[~has_gene]]).sort_index()  # keeps shuffled order
    return capped.head(n)


def build_sample(
    df: pd.DataFrame,
    *,
    cohort: str,
    per_class: int,
    gene_cap: int = config.GENE_CAP,
    seed: int = config.SEED,
    exclude_ids: set[int] | None = None,
    attrition: Attrition | None = None,
) -> pd.DataFrame:
    """Stratified, class-aware sample. `df` must already contain `stratum`."""
    pool = df[df["stratum"].notna()]
    if exclude_ids:
        pool = pool[~pool["variation_id"].isin(exclude_ids)]
    pieces = []
    for s_idx, stratum in enumerate(config.STRATA):
        for label in (1, 0):
            cell = pool[(pool["stratum"] == stratum) & (pool["label"] == label)]
            rng = np.random.default_rng([seed, s_idx, label, 0 if cohort == "primary" else 1])
            pieces.append(sample_cell(cell, per_class, gene_cap, rng))
    out = pd.concat(pieces, ignore_index=True) if pieces else pool.iloc[0:0]
    out = out.assign(cohort=cohort)
    out = out.sort_values(["stratum", "label", "variation_id"]).reset_index(drop=True)
    if attrition is not None:
        attrition.add(f"[{cohort}] final sample", len(out), note=f"per_class<={per_class}, gene_cap={gene_cap}")
    return out[SAMPLE_COLUMNS]


def build_primary_and_onestar(df: pd.DataFrame, attrition: Attrition) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Primary cohort (>=2 stars) and supplementary 1-star-only cohort (design.md deviation 2)."""
    attrition.add("[primary] review status >= 2 stars", int((df["stars"] >= config.PRIMARY_MIN_STARS).sum()))
    primary_pool = df[df["stars"] >= config.PRIMARY_MIN_STARS]
    mapped = primary_pool[primary_pool["stratum"].notna()]
    attrition.add("[primary] mapped to a consequence stratum", len(mapped),
                  note=f"unmapped: {len(primary_pool) - len(mapped)}")
    primary = build_sample(mapped, cohort="primary", per_class=config.PER_CLASS_PER_STRATUM, attrition=attrition)
    if len(primary) > config.MAX_SAMPLE_TOTAL:
        raise AssertionError(f"Primary sample has {len(primary)} variants, above the cap {config.MAX_SAMPLE_TOTAL}")

    onestar_pool = df[df["stars"] == 1]
    onestar = build_sample(
        onestar_pool, cohort="onestar", per_class=config.ONESTAR_PER_CLASS_PER_STRATUM, attrition=attrition,
    )
    return primary, onestar


def pick_dry_run(primary: pd.DataFrame, n: int = config.DRY_RUN_N, seed: int = config.SEED) -> pd.DataFrame:
    """A small, balanced, reproducible subset for the end-to-end dry run."""
    half = n // 2
    rng = np.random.default_rng([seed, 99])
    parts = []
    for label in (1, 0):
        sub = primary[primary["label"] == label].sort_values("variation_id")
        take = min(half, len(sub))
        parts.append(sub.iloc[rng.choice(len(sub), size=take, replace=False)])
    return pd.concat(parts, ignore_index=True)
