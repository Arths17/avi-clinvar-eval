"""Post-hoc diagnostics run after the design.md section 7 sanity gate fired.

EXPLORATORY, NOT PRE-REGISTERED. Output goes to results/diagnostics.json and
results/diagnostics_af_shortcut.csv. Run:  python scripts/diagnostics.py
Needs the cache (scored variants) and data/raw/clinvar_*.vcf.gz.
"""
from __future__ import annotations

import gzip
import json

import numpy as np
import pandas as pd
from sklearn.metrics import roc_auc_score

from avi_eval import analysis, clinvar, config, scoring

OUT = config.RESULTS_DIR
VCF = config.RAW_DIR / f"clinvar_{config.CLINVAR_RELEASE}.vcf.gz"


def auc(d: pd.DataFrame, col: str = "avi_raw") -> float:
    return float("nan") if d["label"].nunique() < 2 else float(roc_auc_score(d["label"], d[col]))


def main() -> None:
    cache = scoring.Cache(config.CACHE_DIR)
    df, _ = analysis.load_scored([config.DATA_DIR / "sample_primary.csv"], cache)
    report: dict = {}

    # A. integrity -----------------------------------------------------------
    report["A_integrity"] = {
        "rows": len(df),
        "variation_id_unique": bool(df["variation_id"].is_unique),
        "coord_keys_unique": bool(df[["chrom", "pos", "ref", "alt"]].drop_duplicates().shape[0] == len(df)),
        "nan_scores": int(df["avi_raw"].isna().sum()),
        "sklearn_pooled_auc": auc(df),
        "sklearn_auc_by_stratum": {s: auc(d) for s, d in df.groupby("stratum")},
    }
    ids = set(df.sample(300, random_state=0)["variation_id"])
    found, af = {}, {}
    all_ids = set(df["variation_id"])
    with gzip.open(VCF, "rt") as fh:
        for line in fh:
            if line[0] == "#":
                continue
            x = line.split("\t", 8)
            if not x[2].isdigit():
                continue
            vid = int(x[2])
            if vid in all_ids:
                d = clinvar._pick_info(x[7], ("CLNSIG", "AF_EXAC", "AF_TGP", "AF_ESP"))
                vals = [float(v) for k, v in d.items() if k.startswith("AF_") and v not in ("", ".")]
                af[vid] = max(vals) if vals else np.nan
                if vid in ids:
                    found[vid] = (clinvar.normalise_chrom(x[0]), int(x[1]), x[3], x[4], d.get("CLNSIG"))
    mism = 0
    for r in df[df["variation_id"].isin(ids)].itertuples():
        c, p, ref, alt, sig = found[r.variation_id]
        mism += int(not ((c, p, ref, alt) == (r.chrom, r.pos, r.ref, r.alt) and clinvar.label_from_clnsig(sig) == r.label))
    report["A_integrity"]["raw_vcf_recheck"] = {"n_checked": len(found), "mismatches": mism}

    # B. label-shuffle null (within stratum) ---------------------------------------
    rng = np.random.default_rng(config.SEED)
    macro = []
    for _ in range(500):
        macro.append(np.mean([roc_auc_score(rng.permutation(d["label"].to_numpy()), d["avi_raw"])
                              for _, d in df.groupby("stratum") if d["label"].nunique() == 2]))
    report["B_shuffle_null_macro_auc"] = {"mean": float(np.mean(macro)), "p2.5": float(np.percentile(macro, 2.5)),
                                          "p97.5": float(np.percentile(macro, 97.5))}

    # C. which feature carries the signal ------------------------------------------
    feat = pd.DataFrame(list(df["features"]), index=df.index)
    top = {}
    for s in config.STRATA:
        m = df["stratum"] == s
        if df.loc[m, "label"].nunique() < 2:
            continue
        a = {c: float(roc_auc_score(df.loc[m, "label"], feat.loc[m, c])) for c in feat.columns if feat.loc[m, c].nunique() > 1}
        top[s] = dict(sorted(a.items(), key=lambda kv: -abs(kv[1] - 0.5))[:4])
    report["C_top_single_feature_auc_by_stratum"] = top

    # D. common-variant ("allele frequency") shortcut ---------------------------------
    df["af_max"] = df["variation_id"].map(af)
    df["has_af"] = df["af_max"].notna()
    rows = []
    for s, d in list(df.groupby("stratum")) + [("POOLED", df)]:
        nf, wf = d[~d["has_af"]], d[d["has_af"]]
        rows.append(
            {
                "stratum": s,
                "share_with_AF_pathogenic": float(d[d["label"] == 1]["has_af"].mean()),
                "share_with_AF_benign": float(d[d["label"] == 0]["has_af"].mean()),
                "auc_all": auc(d), "auc_no_AF": auc(nf), "n_pos_noAF": int((nf["label"] == 1).sum()),
                "n_neg_noAF": int((nf["label"] == 0).sum()), "auc_with_AF": auc(wf),
                "n_pos_AF": int((wf["label"] == 1).sum()), "n_neg_AF": int((wf["label"] == 0).sum()),
            }
        )
    pd.DataFrame(rows).to_csv(OUT / "diagnostics_af_shortcut.csv", index=False)
    report["D_af_shortcut_file"] = "results/diagnostics_af_shortcut.csv"
    report["D_median_af_among_variants_with_af"] = df.groupby("label")["af_max"].median().to_dict()

    # E. leakage via AVI's own training chromosomes -------------------------------------
    ho = df[df["heldout_chrom"]]
    report["E_heldout_chrom_auc"] = {s: auc(d) for s, d in ho.groupby("stratum") if d["label"].nunique() == 2 and len(d) >= 60}

    (OUT / "diagnostics.json").write_text(json.dumps(report, indent=2, default=str))
    print(json.dumps(report, indent=2, default=str)[:3500])


if __name__ == "__main__":
    main()
