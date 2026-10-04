"""Exploratory diagnostic (design.md deviation 7): is AVI's high AUC in the synonymous and
intronic strata concentrated at variants right next to exon boundaries?

For every sampled variant: distance (bp) to the nearest exon edge in GENCODE v46 (basic set).
  0 = the variant is the first/last base of an exon; 1 = the first intronic base next to an exon.
Bands were fixed before looking at any result: 0-2, 3-10, 11-50, 51-200, >200.
AUC is reported only where each class has >= 10 variants (looser than the pre-registered 30,
because this is exploratory). Needs data/raw/gencode.v46.basic.annotation.gtf.gz.
Run:  python scripts/splice_distance.py
"""
from __future__ import annotations

import gzip
from collections import defaultdict

import numpy as np
import pandas as pd
from sklearn.metrics import roc_auc_score

from avi_eval import analysis, config, metrics, scoring

GTF = config.RAW_DIR / "gencode.v46.basic.annotation.gtf.gz"
BANDS = [(0, 2, "0-2"), (3, 10, "3-10"), (11, 50, "11-50"), (51, 200, "51-200"), (201, 10**9, ">200")]
MIN_PER_CLASS = 10


def exon_edges() -> dict[str, np.ndarray]:
    edges: dict[str, list[int]] = defaultdict(list)
    with gzip.open(GTF, "rt") as fh:
        for line in fh:
            if line[0] == "#":
                continue
            f = line.split("\t", 8)
            if f[2] == "exon":
                edges[f[0]].extend((int(f[3]), int(f[4])))  # 1-based inclusive start/end
    return {c: np.unique(np.array(v)) for c, v in edges.items()}


def nearest_edge_distance(edges: np.ndarray, pos: int) -> int:
    i = np.searchsorted(edges, pos)
    cands = []
    if i < len(edges):
        cands.append(abs(int(edges[i]) - pos))
    if i > 0:
        cands.append(abs(pos - int(edges[i - 1])))
    return min(cands)


def band_of(d: int) -> str:
    for lo, hi, name in BANDS:
        if lo <= d <= hi:
            return name
    return ">200"


def auc(d: pd.DataFrame) -> float:
    n1, n0 = int((d["label"] == 1).sum()), int((d["label"] == 0).sum())
    return float(roc_auc_score(d["label"], d["avi_raw"])) if min(n1, n0) >= MIN_PER_CLASS else float("nan")


def main() -> None:
    cache = scoring.Cache(config.CACHE_DIR)
    df, _ = analysis.load_scored([config.DATA_DIR / "sample_primary.csv"], cache)
    edges = exon_edges()
    dist = []
    for r in df.itertuples():
        e = edges.get(f"chr{r.chrom}")
        dist.append(nearest_edge_distance(e, int(r.pos)) if e is not None else -1)
    df["dist_to_exon_edge"] = dist
    print("variants with no exon on their chromosome in the annotation:", int((df["dist_to_exon_edge"] < 0).sum()))
    df = df[df["dist_to_exon_edge"] >= 0].copy()
    df["band"] = df["dist_to_exon_edge"].map(band_of)

    rows = []
    for stratum in ["intronic", "synonymous", "missense", "ALL"]:
        sub = df if stratum == "ALL" else df[df["stratum"] == stratum]
        for _, _, band in [(0, 0, "ALL BANDS")] + BANDS:
            d = sub if band == "ALL BANDS" else sub[sub["band"] == band]
            rows.append(
                {
                    "stratum": stratum, "band": band, "n_pathogenic": int((d["label"] == 1).sum()),
                    "n_benign": int((d["label"] == 0).sum()), "auc_avi": auc(d),
                    "auc_splicing_feature": float("nan") if d.empty else (
                        float(roc_auc_score(d["label"], [f["MERGED_SPLICING"] for f in d["features"]]))
                        if min((d["label"] == 1).sum(), (d["label"] == 0).sum()) >= MIN_PER_CLASS else float("nan")),
                }
            )
    out = pd.DataFrame(rows)
    # gene-clustered bootstrap 95% CI for the AVI AUC of every (stratum, band) cell (same method as the main analysis)
    df = df.reset_index(drop=True)
    groups = {}
    for stratum in ["intronic", "synonymous", "missense", "ALL"]:
        base = np.ones(len(df), bool) if stratum == "ALL" else (df["stratum"] == stratum).to_numpy()
        groups[(stratum, "ALL BANDS")] = base
        for _, _, band in BANDS:
            groups[(stratum, band)] = base & (df["band"] == band).to_numpy()
    res = metrics.bootstrap(df, ["avi_raw"], {f"{a}|{b}": m for (a, b), m in groups.items()}, n_boot=config.N_BOOT)
    ci = {k[1]: metrics.percentile_ci(v) for k, v in res["boot_roc"].items()}
    out["auc_lo"] = [ci[f"{r.stratum}|{r.band}"][0] if r.auc_avi == r.auc_avi else np.nan for r in out.itertuples()]
    out["auc_hi"] = [ci[f"{r.stratum}|{r.band}"][1] if r.auc_avi == r.auc_avi else np.nan for r in out.itertuples()]
    out.to_csv(config.RESULTS_DIR / "diagnostics_splice_distance.csv", index=False)
    # composition: where do pathogenic vs benign variants sit?
    comp = (df[df["stratum"].isin(["intronic", "synonymous"])]
            .groupby(["stratum", "label", "band"]).size().unstack("band", fill_value=0)
            .reindex(columns=[b[2] for b in BANDS], fill_value=0))
    comp.to_csv(config.RESULTS_DIR / "diagnostics_splice_distance_composition.csv")
    pd.set_option("display.width", 200)
    print(out.round(3).to_string(index=False))
    print("\nWhere the sampled variants sit (rows: stratum, label 1=pathogenic):\n", comp.to_string())


if __name__ == "__main__":
    main()
