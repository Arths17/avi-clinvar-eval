"""Rank-based metrics and the gene-clustered bootstrap (design.md section 5).

Nothing here fits a model to any score: only ranks and counts are used.
"""
from __future__ import annotations

from typing import Mapping, Sequence

import numpy as np
import pandas as pd
from scipy.stats import rankdata
from sklearn.metrics import average_precision_score

from . import config


# --------------------------------------------------------------------------- #
# Point estimates
# --------------------------------------------------------------------------- #
def roc_auc(y: np.ndarray, s: np.ndarray) -> float:
    """P(score of a random positive > score of a random negative), ties count 1/2.

    Returns NaN if either class is empty. Equivalent to the Mann-Whitney U statistic / (n1*n0).
    """
    y = np.asarray(y)
    s = np.asarray(s, dtype=float)
    n1 = int((y == 1).sum())
    n0 = int((y == 0).sum())
    if n1 == 0 or n0 == 0:
        return float("nan")
    ranks = rankdata(s)  # average ranks for ties
    return float((ranks[y == 1].sum() - n1 * (n1 + 1) / 2.0) / (n1 * n0))


def avg_precision(y: np.ndarray, s: np.ndarray) -> float:
    """PR-AUC as average precision. NaN if either class is empty."""
    y = np.asarray(y)
    if (y == 1).sum() == 0 or (y == 0).sum() == 0:
        return float("nan")
    return float(average_precision_score(y, s))


def prevalence(y: np.ndarray) -> float:
    y = np.asarray(y)
    return float("nan") if len(y) == 0 else float(np.mean(y == 1))


# --------------------------------------------------------------------------- #
# Baseline B1 and thresholds
# --------------------------------------------------------------------------- #
def severity_baseline(strata: Sequence[str]) -> np.ndarray:
    """Fixed ordinal consequence-severity score (config.SEVERITY). Not fitted to data."""
    return np.array([config.SEVERITY[s] for s in strata], dtype=float)


def threshold_table(df: pd.DataFrame, thresholds: Sequence[int] = config.PHRED_THRESHOLDS) -> pd.DataFrame:
    """Sensitivity / specificity at fixed PHRED cut-offs, pooled and per stratum."""
    rows = []
    groups: list[tuple[str, pd.DataFrame]] = [("pooled", df)] + [(s, g) for s, g in df.groupby("stratum")]
    for name, g in groups:
        pos, neg = g[g["label"] == 1], g[g["label"] == 0]
        for t in thresholds:
            tp = int((pos["avi_phred"] >= t).sum())
            tn = int((neg["avi_phred"] < t).sum())
            rows.append(
                {
                    "group": name, "phred_threshold": t, "n_pos": len(pos), "n_neg": len(neg),
                    "sensitivity": tp / len(pos) if len(pos) else np.nan,
                    "specificity": tn / len(neg) if len(neg) else np.nan,
                }
            )
    return pd.DataFrame(rows)


# --------------------------------------------------------------------------- #
# Gene-clustered bootstrap
# --------------------------------------------------------------------------- #
def cluster_ids(df: pd.DataFrame) -> np.ndarray:
    """Gene symbol as cluster; variants with no gene are their own cluster."""
    return np.array(
        [g if isinstance(g, str) and g else f"__nogene_{v}" for g, v in zip(df["gene"], df["variation_id"])]
    )


def eligible(y: np.ndarray, min_per_class: int = config.MIN_PER_CLASS_FOR_AUC) -> bool:
    y = np.asarray(y)
    return bool((y == 1).sum() >= min_per_class and (y == 0).sum() >= min_per_class)


def bootstrap(
    df: pd.DataFrame,
    scores: Sequence[str],
    groups: Mapping[str, np.ndarray],
    *,
    macro_over: Sequence[str] = (),
    n_boot: int = config.N_BOOT,
    seed: int = config.SEED,
) -> dict:
    """Gene-clustered bootstrap for several scores and several row-subsets at once.

    Parameters
    ----------
    df : rows with `label`, `gene`, `variation_id` and each column in `scores`.
    groups : name -> boolean mask over df rows ("pooled", each stratum, roll-ups, ...).
    macro_over : names in `groups` whose ROC-AUCs are averaged (equal weight) into "macro".

    Returns a dict with point estimates and replicate arrays:
        {"point": {(score, group): (roc, ap)},
         "boot_roc": {(score, group): array}, "boot_ap": {...}, "boot_macro": {score: array}, "point_macro": {score: float}}
    """
    y_all = df["label"].to_numpy()
    cl = cluster_ids(df)
    uniq, inv = np.unique(cl, return_inverse=True)
    rows_by_cluster = [np.flatnonzero(inv == i) for i in range(len(uniq))]
    rng = np.random.default_rng(seed)
    S = {s: df[s].to_numpy(dtype=float) for s in scores}
    masks = {g: np.asarray(m, dtype=bool) for g, m in groups.items()}

    point, boot_roc, boot_ap = {}, {}, {}
    for s in scores:
        for g, m in masks.items():
            point[(s, g)] = (roc_auc(y_all[m], S[s][m]), avg_precision(y_all[m], S[s][m]))
            boot_roc[(s, g)] = np.full(n_boot, np.nan)
            boot_ap[(s, g)] = np.full(n_boot, np.nan)
    boot_macro = {s: np.full(n_boot, np.nan) for s in scores}
    point_macro = {
        s: float(np.nanmean([point[(s, g)][0] for g in macro_over])) if macro_over else float("nan") for s in scores
    }

    n_clusters = len(uniq)
    for b in range(n_boot):
        draw = rng.integers(0, n_clusters, size=n_clusters)
        idx = np.concatenate([rows_by_cluster[i] for i in draw])
        yb = y_all[idx]
        for s in scores:
            sb = S[s][idx]
            macro_vals = []
            for g, m in masks.items():
                mb = m[idx]
                if mb.sum() == 0:
                    continue
                r = roc_auc(yb[mb], sb[mb])
                boot_roc[(s, g)][b] = r
                boot_ap[(s, g)][b] = avg_precision(yb[mb], sb[mb])
                if g in macro_over and not np.isnan(r):
                    macro_vals.append(r)
            if macro_vals:
                boot_macro[s][b] = float(np.mean(macro_vals))
    return {
        "point": point, "boot_roc": boot_roc, "boot_ap": boot_ap,
        "boot_macro": boot_macro, "point_macro": point_macro,
    }


def percentile_ci(values: np.ndarray, level: float = 0.95) -> tuple[float, float]:
    v = values[~np.isnan(values)]
    if len(v) == 0:
        return (float("nan"), float("nan"))
    a = (1 - level) / 2
    return (float(np.quantile(v, a)), float(np.quantile(v, 1 - a)))


def leave_top_genes_out(df: pd.DataFrame, score: str, k: int = 5) -> dict:
    """Point estimates after removing the k genes with the most rows in the sample."""
    counts = df["gene"].dropna().value_counts()
    top = list(counts.index[:k])
    sub = df[~df["gene"].isin(top)]
    return {
        "removed_genes": top,
        "n_remaining": len(sub),
        "roc_auc": roc_auc(sub["label"].to_numpy(), sub[score].to_numpy()),
        "pr_auc": avg_precision(sub["label"].to_numpy(), sub[score].to_numpy()),
    }
