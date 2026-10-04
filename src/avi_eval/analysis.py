"""Turn scored samples into the tables defined in design.md sections 5-7."""
from __future__ import annotations

import json
import logging
import re
from pathlib import Path

import numpy as np
import pandas as pd

from . import config, metrics
from .scoring import Cache, attach_scores

log = logging.getLogger(__name__)

SCORE_AVI = "avi_raw"
SCORE_B1 = "severity_b1"
SCORE_AG = "ag_shap_sum"
SCORE_AM = "alphamissense_shap"


def _norm(name: str) -> str:
    return re.sub(r"[^A-Z0-9]", "", str(name).upper())


def feature_columns(features: dict | None) -> tuple[float, float]:
    """(AlphaMissense SHAP, sum of the ten AlphaGenome-derived SHAPs) from one feature dict.

    Raises if a required feature name cannot be matched - we never guess.
    """
    if features is None:
        return (np.nan, np.nan)
    by_norm = {_norm(k): v for k, v in features.items()}
    try:
        am = by_norm[_norm(config.ALPHAMISSENSE_FEATURE)]
        ag = sum(by_norm[_norm(n)] for n in config.ALPHAGENOME_FEATURES)
    except KeyError as exc:
        raise ValueError(
            f"Feature name {exc} not found in AVI feature breakdown. Names seen: {sorted(features)}. "
            "Check config.ALPHAGENOME_FEATURES against the real output before analysing."
        ) from None
    return (float(am), float(ag))


def load_scored(sample_paths: list[Path], cache: Cache) -> tuple[pd.DataFrame, dict]:
    """Read sample CSVs, join cached scores, add baseline/decomposition columns."""
    sample = pd.concat([pd.read_csv(p) for p in sample_paths], ignore_index=True)
    sample = sample.drop_duplicates("variation_id", keep="first")
    scored = attach_scores(sample, cache)
    n_missing = int(scored[SCORE_AVI].isna().sum())
    info = {"n_sampled": len(scored), "n_missing_score": n_missing}
    scored = scored[scored[SCORE_AVI].notna()].copy()
    scored[SCORE_B1] = metrics.severity_baseline(scored["stratum"])
    parts = [feature_columns(f) for f in scored["features"]]
    scored[SCORE_AM] = [p[0] for p in parts]
    scored[SCORE_AG] = [p[1] for p in parts]
    return scored.reset_index(drop=True), info


def make_groups(df: pd.DataFrame) -> tuple[dict[str, np.ndarray], list[str]]:
    """Row masks: pooled, each stratum, roll-ups, held-out chromosomes. Returns (groups, macro_strata)."""
    groups: dict[str, np.ndarray] = {"pooled": np.ones(len(df), dtype=bool)}
    macro = []
    for s in config.STRATA:
        m = (df["stratum"] == s).to_numpy()
        if m.sum():
            groups[f"stratum:{s}"] = m
            if metrics.eligible(df.loc[m, "label"].to_numpy()):
                macro.append(f"stratum:{s}")
    for name, members in config.ROLLUPS.items():
        groups[f"rollup:{name}"] = df["stratum"].isin(members).to_numpy()
    groups["heldout_chroms:pooled"] = df["heldout_chrom"].to_numpy(dtype=bool)
    for s in config.STRATA:
        m = ((df["stratum"] == s) & df["heldout_chrom"]).to_numpy()
        if m.sum():
            groups[f"heldout_chroms:stratum:{s}"] = m
    return groups, macro


def analyse(df: pd.DataFrame, analysis_set: str, n_boot: int) -> pd.DataFrame:
    """Full metric table for one analysis set (e.g. 'primary' or 'onestar_union')."""
    groups, macro = make_groups(df)
    scores = [SCORE_AVI, SCORE_B1, SCORE_AG, SCORE_AM]
    res = metrics.bootstrap(df, scores, groups, macro_over=macro, n_boot=n_boot)
    rows = []
    for (score, g), (roc, ap) in res["point"].items():
        m = groups[g]
        y = df.loc[m, "label"].to_numpy()
        lo_r, hi_r = metrics.percentile_ci(res["boot_roc"][(score, g)])
        lo_a, hi_a = metrics.percentile_ci(res["boot_ap"][(score, g)])
        rows.append(
            {
                "analysis_set": analysis_set, "score": score, "group": g,
                "n_pos": int((y == 1).sum()), "n_neg": int((y == 0).sum()),
                "prevalence": metrics.prevalence(y),
                "eligible_for_auc": metrics.eligible(y),
                "roc_auc": roc, "roc_lo": lo_r, "roc_hi": hi_r,
                "pr_auc": ap, "pr_lo": lo_a, "pr_hi": hi_a,
            }
        )
    for score in scores:
        lo, hi = metrics.percentile_ci(res["boot_macro"][score])
        rows.append(
            {"analysis_set": analysis_set, "score": score, "group": "macro_eligible_strata",
             "n_pos": np.nan, "n_neg": np.nan, "prevalence": np.nan, "eligible_for_auc": True,
             "roc_auc": res["point_macro"][score], "roc_lo": lo, "roc_hi": hi,
             "pr_auc": np.nan, "pr_lo": np.nan, "pr_hi": np.nan}
        )
    # paired difference AVI - B1 on the pooled set (same bootstrap draws)
    for g in ("pooled", "heldout_chroms:pooled"):
        d = res["boot_roc"][(SCORE_AVI, g)] - res["boot_roc"][(SCORE_B1, g)]
        lo, hi = metrics.percentile_ci(d)
        rows.append(
            {"analysis_set": analysis_set, "score": "avi_minus_b1", "group": g,
             "n_pos": np.nan, "n_neg": np.nan, "prevalence": np.nan, "eligible_for_auc": True,
             "roc_auc": res["point"][(SCORE_AVI, g)][0] - res["point"][(SCORE_B1, g)][0],
             "roc_lo": lo, "roc_hi": hi, "pr_auc": np.nan, "pr_lo": np.nan, "pr_hi": np.nan}
        )
    return pd.DataFrame(rows)


def sanity_flags(table: pd.DataFrame) -> dict:
    """Design.md section 7 gate: flag implausible results BEFORE any interpretation."""
    t = table[(table["analysis_set"] == "primary") & (table["score"] == SCORE_AVI)]
    flags = []
    pooled = t[t["group"] == "pooled"].iloc[0]
    if pooled["roc_auc"] > config.SANITY_AUC_HIGH or pooled["roc_auc"] < config.SANITY_AUC_LOW:
        flags.append(f"pooled ROC-AUC {pooled['roc_auc']:.3f} is outside [{config.SANITY_AUC_LOW}, {config.SANITY_AUC_HIGH}]")
    for _, r in t[t["group"].str.startswith("stratum:") & t["eligible_for_auc"]].iterrows():
        if r["roc_auc"] > config.SANITY_AUC_HIGH or r["roc_auc"] < 0.5 - 0.02:
            flags.append(f"{r['group']} ROC-AUC {r['roc_auc']:.3f} looks implausible")
    return {"needs_investigation": bool(flags), "flags": flags}


def run_all(sample_primary: Path, sample_onestar: Path | None, cache: Cache, out_dir: Path, n_boot: int) -> dict:
    out_dir.mkdir(parents=True, exist_ok=True)
    primary, info = load_scored([sample_primary], cache)
    log.info("primary: %s", info)
    tables = [analyse(primary, "primary", n_boot)]

    union = None
    if sample_onestar is not None and sample_onestar.exists():
        union, info_u = load_scored([sample_primary, sample_onestar], cache)
        log.info("onestar_union: %s", info_u)
        tables.append(analyse(union, "onestar_union", n_boot))
    table = pd.concat(tables, ignore_index=True)
    table.to_csv(out_dir / "metrics.csv", index=False)

    counts = (
        primary.groupby(["stratum", "label"]).size().unstack(fill_value=0)
        .rename(columns={0: "n_benign", 1: "n_pathogenic"}).reindex(list(config.STRATA)).fillna(0).astype(int)
    )
    counts.to_csv(out_dir / "class_counts_primary.csv")
    metrics.threshold_table(primary).to_csv(out_dir / "thresholds_primary.csv", index=False)
    ltgo = metrics.leave_top_genes_out(primary, SCORE_AVI)
    ties = {
        "n_unique_avi_raw": int(primary[SCORE_AVI].nunique()),
        "n_rows": len(primary),
        "n_unique_phred": int(primary["avi_phred"].nunique()),
    }
    sanity = sanity_flags(table)
    summary = {"primary_load": info, "leave_top5_genes_out": ltgo, "ties": ties, "sanity": sanity,
               "n_boot": n_boot, "seed": config.SEED}
    (out_dir / "summary.json").write_text(json.dumps(summary, indent=2, default=str))
    return {"table": table, "primary": primary, "union": union, "counts": counts, "summary": summary}
