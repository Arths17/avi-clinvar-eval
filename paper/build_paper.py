"""Build paper/preprint.md (and HTML/PDF) from the real result files.

Every number in the paper is read from results/*.csv|json at build time; none is typed by hand.
Run from the repository root:   python paper/build_paper.py
PDF export uses Google Chrome in headless mode if it is installed.
"""
from __future__ import annotations

import json
import shutil
import subprocess
from pathlib import Path

import markdown
import pandas as pd

R = Path("results")
P = Path("paper")
met = pd.read_csv(R / "metrics.csv")
att = pd.read_csv(R / "attrition.csv").set_index("step")
summ = json.loads((R / "summary.json").read_text())
diag = json.loads((R / "diagnostics.json").read_text())
afd = pd.read_csv(R / "diagnostics_af_shortcut.csv")
spl = pd.read_csv(R / "diagnostics_splice_distance.csv")
thr = pd.read_csv(R / "thresholds_primary.csv")
man = json.loads((R / "run_manifest.json").read_text())
counts = pd.read_csv(R / "class_counts_primary.csv").set_index("stratum")

LABEL = {
    "nonsense_start_stop": "Nonsense / start-lost / stop-lost", "splice_canonical": "Canonical splice site (±1–2 bp)",
    "missense": "Missense", "synonymous": "Synonymous", "splice_region": "Splice region (empty, see Methods)",
    "intronic": "Intronic", "utr5": "5′ UTR", "utr3": "3′ UTR", "other_noncoding": "Other non-coding",
}
ORDER = ["missense", "synonymous", "nonsense_start_stop", "splice_canonical", "intronic", "utr5", "utr3", "other_noncoding"]


def row(aset: str, score: str, group: str) -> pd.Series:
    m = met[(met.analysis_set == aset) & (met.score == score) & (met.group == group)]
    assert len(m) == 1, (aset, score, group, len(m))
    return m.iloc[0]


def ci(r: pd.Series, a="roc", d=3) -> str:
    return f"{r[a + '_auc']:.{d}f} ({r[a + '_lo']:.{d}f}–{r[a + '_hi']:.{d}f})" if a == "roc" else \
        f"{r['pr_auc']:.{d}f} ({r['pr_lo']:.{d}f}–{r['pr_hi']:.{d}f})"


def md_table(df: pd.DataFrame) -> str:
    head = "| " + " | ".join(df.columns) + " |\n|" + "|".join("---" for _ in df.columns) + "|\n"
    return head + "\n".join("| " + " | ".join(str(v) for v in r) + " |" for r in df.itertuples(index=False)) + "\n"


A = lambda g: row("primary", "avi_raw", g)       # noqa: E731
B1 = lambda g: row("primary", "severity_b1", g)  # noqa: E731
pooled, macro = A("pooled"), A("macro_eligible_strata")
delta = row("primary", "avi_minus_b1", "pooled")
n_primary = int(att.loc["[primary] final sample", "n_remaining"])
n_one = int(att.loc["[onestar] final sample", "n_remaining"])
n_pos, n_neg = int(pooled.n_pos), int(pooled.n_neg)
heldout = A("heldout_chroms:pooled")
ltgo = summ["leave_top5_genes_out"]

# ---- tables ---------------------------------------------------------------------------
t_attr = md_table(pd.DataFrame(
    [{"Step": s, "Variants remaining": f"{int(r.n_remaining):,}"} for s, r in att.iterrows()]))

rows = []
for s in ORDER:
    g = f"stratum:{s}"
    r = A(g)
    elig = bool(r.eligible_for_auc)
    rows.append({
        "Category": LABEL[s], "Pathogenic / benign": f"{int(r.n_pos)} / {int(r.n_neg)}",
        "ROC-AUC (95% CI)": ci(r) if elig else "not computed (<30 in a class)",
        "PR-AUC (95% CI)": ci(r, "pr") if elig else "—",
        "Share pathogenic": f"{r.prevalence:.2f}",
    })
rows.append({"Category": "**Pooled**", "Pathogenic / benign": f"{n_pos} / {n_neg}",
             "ROC-AUC (95% CI)": ci(pooled), "PR-AUC (95% CI)": ci(pooled, "pr"), "Share pathogenic": f"{pooled.prevalence:.2f}"})
rows.append({"Category": "**Macro-average of eligible categories**", "Pathogenic / benign": "—",
             "ROC-AUC (95% CI)": ci(macro), "PR-AUC (95% CI)": "—", "Share pathogenic": "—"})
t_main = md_table(pd.DataFrame(rows))

rows = []
for g, name in [("rollup:coding", "Coding (nonsense/start/stop, missense, synonymous)"),
                ("rollup:splice", "Splice (canonical sites only)"), ("rollup:non_coding", "Non-coding (intronic, UTR, other)")]:
    r = A(g)
    rows.append({"Roll-up": name, "Pathogenic / benign": f"{int(r.n_pos)} / {int(r.n_neg)}", "ROC-AUC (95% CI)": ci(r)})
t_roll = md_table(pd.DataFrame(rows))

rows = []
for s in ORDER:
    ra, rb = row("primary", "alphamissense_shap", f"stratum:{s}"), row("primary", "ag_shap_sum", f"stratum:{s}")
    if not bool(A(f"stratum:{s}").eligible_for_auc):
        continue
    rows.append({"Category": LABEL[s], "AVI": f"{A(f'stratum:{s}').roc_auc:.3f}",
                 "AlphaMissense contribution only": ci(ra), "AlphaGenome-derived contributions only": ci(rb)})
t_dec = md_table(pd.DataFrame(rows))

rows = []
for s in ORDER:
    g = f"heldout_chroms:stratum:{s}"
    if g in set(met.group):
        r = A(g)
        if bool(r.eligible_for_auc):
            rows.append({"Category": LABEL[s], "Pathogenic / benign": f"{int(r.n_pos)} / {int(r.n_neg)}",
                         "ROC-AUC (95% CI)": ci(r), "All chromosomes": f"{A('stratum:' + s).roc_auc:.3f}"})
t_held = md_table(pd.DataFrame(rows))

rows = []
for s in ORDER:
    g = f"stratum:{s}"
    pr, un = A(g), row("onestar_union", "avi_raw", g)
    if bool(un.eligible_for_auc):
        rows.append({"Category": LABEL[s], "Pathogenic / benign (≥1 star)": f"{int(un.n_pos)} / {int(un.n_neg)}",
                     "ROC-AUC ≥1 star (95% CI)": ci(un), "ROC-AUC ≥2 stars (primary)": f"{pr.roc_auc:.3f}"})
un_pool, un_macro = row("onestar_union", "avi_raw", "pooled"), row("onestar_union", "avi_raw", "macro_eligible_strata")
rows.append({"Category": "**Pooled**", "Pathogenic / benign (≥1 star)": f"{int(un_pool.n_pos)} / {int(un_pool.n_neg)}",
             "ROC-AUC ≥1 star (95% CI)": ci(un_pool), "ROC-AUC ≥2 stars (primary)": f"{pooled.roc_auc:.3f}"})
t_sens = md_table(pd.DataFrame(rows))

tp = thr[thr.group == "pooled"]
t_thr = md_table(pd.DataFrame([{"AVI PHRED threshold": f"≥ {int(r.phred_threshold)}", "Sensitivity (pathogenic called)": f"{r.sensitivity:.3f}",
                                "Specificity (benign not called)": f"{r.specificity:.3f}"} for r in tp.itertuples()]))

rows = []
for r in afd.itertuples():
    rows.append({"Category": LABEL.get(r.stratum, "**Pooled**"),
                 "Share with a frequency: pathogenic / benign": f"{r.share_with_AF_pathogenic:.2f} / {r.share_with_AF_benign:.2f}",
                 "AUC, no frequency (path./benign n)": f"{r.auc_no_AF:.3f} ({r.n_pos_noAF}/{r.n_neg_noAF})",
                 "AUC, has frequency (path./benign n)": f"{r.auc_with_AF:.3f} ({r.n_pos_AF}/{r.n_neg_AF})"})
t_af = md_table(pd.DataFrame(rows))

rows = []
for r in spl[spl.stratum.isin(["intronic", "synonymous"]) & (spl.band != "ALL BANDS") & ((spl.n_pathogenic + spl.n_benign) > 0)].itertuples():
    auc_txt = f"{r.auc_avi:.3f} ({r.auc_lo:.3f}–{r.auc_hi:.3f})" if r.auc_avi == r.auc_avi else "not computed (<10 in a class)"
    rows.append({"Category": LABEL[r.stratum], "Distance to nearest exon edge (bp)": r.band,
                 "Pathogenic / benign": f"{r.n_pathogenic} / {r.n_benign}", "AVI ROC-AUC (95% CI)": auc_txt})
t_spl = md_table(pd.DataFrame(rows))

sm = lambda s, lab: int(counts.loc[s, lab])  # noqa: E731

from avi_eval import config as _cfg, scoring as _sc  # noqa: E402
_cache = _sc.Cache(_cfg.CACHE_DIR)
_recs = list(_cache.records().values())
n_cached = len(_recs)
n_errors = sum(1 for _ in open(_cache.error_path)) if _cache.error_path.exists() else 0
_off = [r["avi_raw"] - sum(r["features"].values()) for r in _recs]
off_mean, off_min, off_max = sum(_off) / len(_off), min(_off), max(_off)
flag = "; ".join(summ["sanity"]["flags"]) if summ["sanity"]["flags"] else "none"

# ---- the paper ----------------------------------------------------------------------------
TXT = f"""
# How well does the AlphaGenome Variant Impact (AVI) score separate pathogenic from benign human genetic variants? An independent evaluation on ClinVar

**Atharv Ranjan** — Independent researcher, no institutional affiliation  
Preprint, 4 October 2026. Not peer reviewed.  
Cite as: Ranjan A. (2026). Zenodo. doi:10.5281/zenodo.23134413

> **Research use only. This work makes no clinical or diagnostic claims and is not medical advice.** The AVI score and the other AlphaGenome outputs analysed here are for theoretical modelling; they are not intended, validated or approved for clinical use.

## Abstract

**Background.** The AlphaGenome Atlas provides a precomputed AlphaGenome Variant Impact (AVI) score for every possible single-nucleotide variant (SNV) in the human genome. Its developers report strong performance for separating pathogenic from benign ClinVar variants. **Aim.** We evaluated AVI independently on a reproducible, pre-specified sample of ClinVar SNVs. **Methods.** From the ClinVar GRCh38 VCF (file dated 28 September 2026) we sampled {n_primary:,} SNVs with ≥2-star review status ({n_pos:,} pathogenic or likely pathogenic, {n_neg:,} benign or likely benign), stratified by molecular consequence and capped at 10 variants per gene per class and category, retrieved their AVI scores from the Atlas API, and computed ROC-AUC and PR-AUC with gene-clustered bootstrap 95% confidence intervals. **Results.** In this sample, pooled ROC-AUC was {ci(pooled)} and PR-AUC {ci(pooled, 'pr')}; a fixed consequence-type ranking alone reached {B1('pooled').roc_auc:.3f}. The macro-average over categories with at least 30 variants of each class was {ci(macro)}. ROC-AUC was highest for missense ({A('stratum:missense').roc_auc:.3f}), synonymous ({A('stratum:synonymous').roc_auc:.3f}) and intronic ({A('stratum:intronic').roc_auc:.3f}) variants, and lower for canonical splice sites ({A('stratum:splice_canonical').roc_auc:.3f}), nonsense/start/stop variants ({A('stratum:nonsense_start_stop').roc_auc:.3f}) and 5′ UTR variants ({A('stratum:utr5').roc_auc:.3f}). In the synonymous and intronic categories, most pathogenic variants lay close to exon boundaries, so those categories are probably easier than a random sample of such variants would be. **Conclusions.** AVI separated the two ClinVar classes well in this sample, but ClinVar-specific confounds and possible information overlap with ClinVar (which we cannot exclude) limit what can be concluded; the results are not evidence of clinical utility.

## 1. Introduction

Most human genetic variation lies outside protein-coding sequence (about 98% according to [2]), and predicting which variants disrupt function is a central open problem in genetics. AlphaGenome is a deep-learning model from Google DeepMind that predicts many functional genomic signals from DNA sequence and can score the effect of a variant on them [1]. The AlphaGenome Atlas [2] precomputes these scores for roughly 9 billion possible SNVs in the human reference genome (GRCh38) and combines them with protein-level and conservation features into a single supervised score, the **AlphaGenome Variant Impact (AVI) score**. According to the Atlas preprint, AVI was trained to separate observed human variants with a filtering allele frequency above versus below 0.1% (proxy neutral versus proxy impactful variants), not to predict ClinVar labels, and was then benchmarked, among other datasets, on ClinVar [2].

ClinVar is NCBI's public archive of clinical-laboratory and researcher assertions about the significance of variants [3]. Distinguishing pathogenic from benign ClinVar variants is a common benchmark for variant-effect predictors, but it is also easy to over-interpret: ClinVar labels are not experimental ground truth, they are unevenly distributed across genes and consequence types, and a predictor may have seen related information during development.

Independent evaluation matters because the developers of a score have both the most knowledge and the strongest incentive to present it favourably, and because choices such as which review-status filter, which variant types and which sampling to use can change reported performance. Here a sampling and analysis plan was written and frozen before any score was retrieved (`design.md` in the code repository), and every later change is logged. This paper reports one independent researcher's replication on a fresh sample, with an emphasis on how performance varies by variant category and on what could inflate it.

*Related work.* We rely on the sources listed in the References. The AlphaGenome paper [1], the Atlas preprint [2], the terms, the API documentation and the ClinVar and GENCODE documentation pages were read directly; the ClinVar article in [3] is cited as it appears in the reference list of [1] and was not itself opened. The Atlas preprint [2] reports AVI performance on ClinVar as area under the precision-recall curve (AUPRC) at ClinVar's natural class balance, stratified by consequence category. Our sample is deliberately close to class-balanced within categories, so its PR-AUC values are **not comparable** to those figures and we make no direct comparison with them. We did not evaluate any other predictor.

## 2. Data and methods

### 2.1 Pre-registration and deviations
The design (variant source, label rules, sample size, metrics, baselines, confound checks) was fixed in `design.md` before any AVI score was retrieved. Seven dated deviations were then logged; numbers 1–5 were recorded before any AVI score was retrieved and numbers 6–7 after the results had been seen: (1) the ClinVar file was `clinvar_20260928.vcf.gz`, because ClinVar's GRCh38 VCFs are dated roughly weekly and no 1 October file existed; (2) the ≥1-star sensitivity analysis uses the primary sample plus a small supplementary sample of exactly-1-star variants, to stay under the API-call cap; (3) the gene used for capping and clustering is the first symbol in ClinVar's `GENEINFO` field; (4) ClinVar's consequence field contains no splice-region term, so that planned category is empty, and variants with no consequence annotation were excluded; (5) three non-coding categories contain few pathogenic variants (see Results); (6) after the pre-specified sanity gate flagged AUC > 0.97 in two categories, exploratory diagnostics were run; (7) one of those diagnostics used GENCODE exon annotations. Deviations 6 and 7 are **exploratory and not pre-registered**.

### 2.2 Variants and labels
We used the ClinVar GRCh38 VCF `clinvar_20260928.vcf.gz` (MD5 `{man['clinvar_md5']}`, downloaded from NCBI). Only SNVs on chromosomes 1–22 and X were kept. A variant was **positive** if its aggregate germline classification was exactly Pathogenic, Likely pathogenic or Pathogenic/Likely pathogenic, and **negative** if exactly Benign, Likely benign or Benign/Likely benign. Variants of uncertain significance, conflicting classifications, risk factor, drug response, low-penetrance terms, or any record mixing classes were excluded. ClinVar's *review status* (a 0–4 star rating of how well-supported a classification is) was required to be ≥2 stars (several submitters in agreement, an expert panel, or a practice guideline) for the primary analysis, because single-submitter assertions are the noisiest.

{t_attr}
*Table 1. Number of variants remaining after each filtering step (from `results/attrition.csv`).*

### 2.3 Consequence categories and sampling
Each variant was assigned a molecular-consequence category from the Sequence Ontology terms in ClinVar's `MC` field, using a fixed priority order: nonsense/start-lost/stop-lost; canonical splice donor/acceptor; missense; synonymous; intronic; 5′ UTR; 3′ UTR; other non-coding. (Frameshift and other non-SNV types do not occur here. ClinVar's `MC` field contained no splice-region term, so that planned category is empty; {int(att.loc['[primary] review status >= 2 stars', 'n_remaining']) - int(att.loc['[primary] mapped to a consequence stratum', 'n_remaining']):,} otherwise eligible variants with no consequence annotation were excluded.) For each category and class we drew a random sample (fixed seed {summ['seed']}) of up to 250 variants after capping at 10 variants per gene, giving {n_primary:,} primary variants. A supplementary sample of up to 60 variants per class and category of exactly-1-star variants ({n_one:,} variants) was used only for a sensitivity analysis. Each variant was also flagged by whether it lies on a chromosome that AVI used neither for training nor for validation (3, 6, 9, 12, 16, 18, 19, 21; from the Atlas preprint's Methods).

### 2.4 Scores
For each sampled variant we requested the scorers `AVI_SCORE` and `AVI_SCORE_FEATURE_IMPORTANCE` from the AlphaGenome Atlas API (Python package `alphagenome` {man['alphagenome']}) using GRCh38 coordinates and 1-based positions. The primary score is the raw AVI score; higher means more predicted impact. The API also returns a quantile, from which we computed a PHRED-scaled value (−10·log10(1 − quantile); not checked against the Atlas website in this study). The feature breakdown gives each of 18 input features' additive contribution to the score. We grouped these into the AlphaMissense contribution and the sum of the ten AlphaGenome-derived contributions. Across all {n_cached:,} scored variants, the raw score minus the sum of the 18 contributions averaged {off_mean:.3f} (range {off_min:.3f} to {off_max:.3f}), so the contributions do not sum exactly to the raw score and we treat them as ranking scores, not an exact decomposition. {n_cached:,} of {n_primary + n_one:,} requested variants were scored and {n_errors} errors were logged; responses were cached on disk, requests were limited to 4 concurrent workers with backoff, and a hard cap of 6,000 attempts was enforced ({man.get('api_calls_used', 'n/a')} used).

### 2.5 Metrics and statistics
**ROC-AUC** is the probability that a randomly chosen pathogenic variant receives a higher score than a randomly chosen benign one (0.5 is chance, 1.0 is perfect). **PR-AUC** (average precision) measures how enriched the top-ranked variants are for pathogenic ones; its chance level equals the share of pathogenic variants, which we report next to every value. Confidence intervals are 95% percentile intervals from a bootstrap of {summ['n_boot']:,} replicates that **resamples genes, not variants**, because variants within a gene are not independent. A category receives an AUC only if it has at least 30 variants of each class; otherwise it is described but not scored. The macro-average is the unweighted mean of the eligible categories' ROC-AUCs. No model was trained or fitted on any AlphaGenome output; only rank-based metrics and counts were computed.

### 2.6 Baseline and secondary analyses
*Baseline B1* is a fixed ordinal ranking by consequence type (nonsense/start/stop 5 > canonical splice 4 > missense 3 > synonymous 1 > all non-coding 0; fixed before scoring and never tuned). *Decomposition*: AUC of the AlphaMissense contribution alone and of the AlphaGenome-derived contributions alone. *Held-out chromosomes*: the primary analysis restricted to chromosomes AVI did not train or validate on. *Sensitivity*: primary sample plus 1-star variants. *Fixed thresholds*: sensitivity and specificity at PHRED 10, 20 and 30. *Leave-top-genes-out*: the five most frequent genes removed. We did not evaluate other predictors, so **no comparison with any other method is made**.

### 2.7 Exploratory diagnostics (not pre-registered)
After the sanity gate (pooled AUC outside [0.55, 0.97] or any category AUC above 0.97) flagged two categories, we (a) re-derived AUCs with scikit-learn and re-read 300 sampled variants from the raw ClinVar file, (b) ran a within-category label-shuffle null, (c) computed single-feature AUCs, (d) split variants by whether ClinVar lists a population allele frequency (ExAC, 1000 Genomes or ESP), and (e) computed each variant's distance to the nearest exon edge using GENCODE release 46 (basic set, GRCh38) in the bands 0–2, 3–10, 11–50, 51–200 and >200 bp, with AUC reported where each class had at least 10 variants.

### 2.8 Reproducibility
Code, design, tests and the list of sampled ClinVar variants are in the repository (Section 8). Python {man['python']}, pandas {man['pandas']}, numpy {man['numpy']}, scikit-learn {man['scikit-learn']}; random seed {summ['seed']}.

## 3. Results

### 3.1 Overall and by category
{t_main}
*Table 2. AVI discrimination of pathogenic from benign ClinVar SNVs, primary sample (≥2 stars). CIs are gene-clustered bootstrap 95% intervals. "Share pathogenic" is the chance level of the PR-AUC.*

In this sample, AVI achieved a pooled ROC-AUC of {ci(pooled)} and PR-AUC of {ci(pooled, 'pr')} (chance {pooled.prevalence:.2f}). The pooled value is partly produced by the mix of categories, because some categories are almost entirely one class in ClinVar (for example, most nonsense variants are pathogenic and most synonymous variants are benign). The fixed consequence-type baseline alone reached a pooled ROC-AUC of {ci(B1('pooled'))}; AVI exceeded it by {delta.roc_auc:.3f} (95% CI {delta.roc_lo:.3f}–{delta.roc_hi:.3f}). Within a category the baseline is constant by construction, so the more informative summary is the macro-average of {ci(macro)}.

![Pooled ROC curves for AVI and the consequence baseline](figures/roc_pooled.png)
*Figure 1. Pooled ROC curves, primary sample.*

![AVI ROC-AUC by consequence category](figures/auc_by_stratum.png)
*Figure 2. ROC-AUC by category with gene-clustered 95% intervals (categories with ≥30 variants per class).*

{t_roll}
*Table 3. Roll-ups (exploratory grouping of the categories above).*

Two categories have too few pathogenic variants for an AUC: 3′ UTR ({int(A('stratum:utr3').n_pos)} pathogenic) and other non-coding ({int(A('stratum:other_noncoding').n_pos)}). A third, 5′ UTR, has {int(A('stratum:utr5').n_pos)} pathogenic variants and a wide interval (ROC-AUC {ci(A('stratum:utr5'))}). Performance on non-coding variants other than intronic therefore cannot be characterised from this sample.

![Class balance by category](figures/class_counts.png)
*Figure 3. Number of pathogenic and benign variants in each category.*

### 3.2 Which part of AVI carries the signal
{t_dec}
*Table 4. ROC-AUC of AVI and of its two groups of contributions (eligible categories). The "contribution only" columns are not stand-alone predictors.*

For missense variants, the AlphaMissense contribution alone ({ci(row('primary', 'alphamissense_shap', 'stratum:missense'))}) performed about as well as the full score ({A('stratum:missense').roc_auc:.3f}). For synonymous and intronic variants the AlphaMissense contribution carried no signal (~0.5) and the AlphaGenome-derived contributions carried nearly all of it.

### 3.3 Held-out chromosomes, sensitivity and thresholds
{t_held}
*Table 5. Primary analysis restricted to chromosomes AVI did not use for training or validation.*

On those chromosomes the pooled ROC-AUC was {ci(heldout)} ({int(heldout.n_pos)} pathogenic, {int(heldout.n_neg)} benign), close to the all-chromosome value. Removing the five most frequent genes ({', '.join(ltgo['removed_genes'])}) left a pooled ROC-AUC of {ltgo['roc_auc']:.3f}.

{t_sens}
*Table 6. Sensitivity analysis including 1-star variants (primary sample plus the supplementary sample; different class and category mix).*

With 1-star variants added, the pooled ROC-AUC was {ci(un_pool)} and the macro-average {ci(un_macro)}.

{t_thr}
*Table 7. Pooled sensitivity and specificity at fixed AVI PHRED thresholds (not tuned; this sample's class mix is not representative of any clinical population).*

### 3.4 Sanity gate and exploratory diagnostics
The pre-specified sanity gate flagged: {flag}. Investigation found no bug: our AUCs matched scikit-learn, 300 re-read variants showed {diag['A_integrity']['raw_vcf_recheck']['mismatches']} coordinate or label mismatches, and shuffling labels within categories gave a macro-AUC of {diag['B_shuffle_null_macro_auc']['mean']:.3f} (95% range {diag['B_shuffle_null_macro_auc']['p2.5']:.3f}–{diag['B_shuffle_null_macro_auc']['p97.5']:.3f}). In the synonymous and intronic categories the AUC of the splicing feature alone was {diag['C_top_single_feature_auc_by_stratum']['synonymous']['MERGED_SPLICING']:.3f} and {diag['C_top_single_feature_auc_by_stratum']['intronic']['MERGED_SPLICING']:.3f}.

**Allele-frequency shortcut.** ClinVar benign variants far more often carry a population allele frequency than pathogenic ones, and AVI was trained on a common-versus-rare contrast, so we checked whether this explains performance.

{t_af}
*Table 8. Share of variants with a listed ExAC/1000 Genomes/ESP allele frequency, and AUC with and without one (exploratory).*

AUC remained high among variants with no listed frequency (for example intronic {afd.loc[afd.stratum == 'intronic', 'auc_no_AF'].iloc[0]:.3f}, synonymous {afd.loc[afd.stratum == 'synonymous', 'auc_no_AF'].iloc[0]:.3f}), so the frequency shortcut does not explain the high values, although the benign groups without a frequency are small. For nonsense and canonical-splice variants, AUC was much higher without a listed frequency ({afd.loc[afd.stratum == 'nonsense_start_stop', 'auc_no_AF'].iloc[0]:.3f} and {afd.loc[afd.stratum == 'splice_canonical', 'auc_no_AF'].iloc[0]:.3f}) than with one ({afd.loc[afd.stratum == 'nonsense_start_stop', 'auc_with_AF'].iloc[0]:.3f} and {afd.loc[afd.stratum == 'splice_canonical', 'auc_with_AF'].iloc[0]:.3f}); we did not investigate why.

**Proximity to exon boundaries.** Among synonymous variants, {int(spl[(spl.stratum == 'synonymous') & (spl.band == '0-2')].n_pathogenic.iloc[0])} of {sm('synonymous', 'n_pathogenic')} pathogenic variants lay within 0–2 bp of an exon edge, against {int(spl[(spl.stratum == 'synonymous') & (spl.band == '0-2')].n_benign.iloc[0])} of {sm('synonymous', 'n_benign')} benign ones; among intronic variants, {int(spl[(spl.stratum == 'intronic') & (spl.band == '3-10')].n_pathogenic.iloc[0])} of {sm('intronic', 'n_pathogenic')} pathogenic variants lay 3–10 bp from an exon, against {int(spl[(spl.stratum == 'intronic') & (spl.band == '3-10')].n_benign.iloc[0])} of {sm('intronic', 'n_benign')} benign ones. In these categories the ClinVar label is therefore strongly associated with splice-site proximity, a setting that may favour a splicing-based predictor (not tested here).

{t_spl}
*Table 9. AVI ROC-AUC within distance-to-exon bands (exploratory; intronic variants at 1–2 bp are canonical splice sites and form a separate category).*

AVI discriminated within distance bands as well (for example intronic 11–50 bp {ci(pd.Series({'roc_auc': spl[(spl.stratum == 'intronic') & (spl.band == '11-50')].auc_avi.iloc[0], 'roc_lo': spl[(spl.stratum == 'intronic') & (spl.band == '11-50')].auc_lo.iloc[0], 'roc_hi': spl[(spl.stratum == 'intronic') & (spl.band == '11-50')].auc_hi.iloc[0]}))}), but several bands contain few variants of one class and their intervals are wide.

## 4. Discussion

In this sample, AVI ranked ClinVar pathogenic variants above benign ones well in most categories, with ROC-AUC between {A('stratum:utr5').roc_auc:.2f} (5′ UTR, wide interval) and {A('stratum:intronic').roc_auc:.2f} (intronic). Several observations temper this. First, the pooled value is largely a function of the category mix, as the consequence baseline shows. Second, the highest values occur in categories where ClinVar's pathogenic variants are concentrated near exon boundaries, and the diagnostics show the labels there are associated with splice-site proximity. Third, in missense variants the score was matched by its AlphaMissense component, so the sequence model contributed little additional discrimination in that category in this sample (we did not test this formally). Fourth, performance was lower for canonical splice-site and nonsense/start/stop variants; for both, AUC was higher among variants without a listed population frequency than among those with one, and we did not investigate why. Finally, ClinVar's sparse non-coding pathogenic variants outside introns made a meaningful assessment of UTR and other non-coding categories impossible here. We make no claim about how AVI compares with any other predictor.

## 5. Limitations

* **Possible overlap with ClinVar (circularity).** AVI was not trained on ClinVar labels [2], but its inputs include loss-of-function flags derived partly from ClinVar descriptors and AlphaMissense; the developers also evaluated AVI on ClinVar. We did not verify whether ClinVar was used to train or calibrate AlphaMissense or the underlying AlphaGenome model, and the AlphaGenome paper lists ClinVar among its data sources [1]. Our restriction to chromosomes unused by AVI addresses only AVI's own training split. **We cannot exclude information overlap with ClinVar, so these results are not a fully independent test.**
* **Label quality and ClinVar bias.** ClinVar labels are assertions, not experimental truth; they over-represent well-studied genes and variant types, and may over-represent some ancestry groups. The ≥2-star filter may increase this bias while reducing label noise. Benign classification often relies partly on population frequency; in our sample benign variants far more often had a listed frequency than pathogenic ones (Table 8).
* **Sampling.** Class-balanced, gene-capped, stratified sampling means results describe this sample, not ClinVar or the genome. The pooled values depend on our category quotas. Per-gene results are not reported.
* **Small groups.** Several categories (3′ UTR, other non-coding, 5′ UTR) have very few pathogenic variants; many subgroup analyses have wide intervals.
* **Exploratory analyses** (Section 3.4) were chosen after seeing results and carry no pre-registered status. The distance-to-exon analysis uses the nearest exon edge in a transcript set without regard to gene.
* **Single release, single build, SNVs only.** One ClinVar release, GRCh38, no indels, no other predictors or baselines beyond a trivial consequence ranking.
* **Score details.** The PHRED convention was not verified against the Atlas website; the primary analysis uses the raw score and is unaffected. The feature contributions do not sum exactly to the raw score.
* **Research use only.** Nothing here supports clinical interpretation of any variant.

## 6. Conclusion

In a reproducible, pre-specified sample of {n_primary:,} ClinVar SNVs, the AVI score achieved a pooled ROC-AUC of {pooled.roc_auc:.2f} and a category-averaged ROC-AUC of {macro.roc_auc:.2f}, well above a trivial consequence-type ranking, with the highest values in missense, synonymous and intronic variants. Much of the difficulty and the apparent ease of different categories is tied to how ClinVar labels are distributed, and independence from ClinVar cannot be guaranteed. These findings describe one sample of one database and should not be taken to show clinical utility.

## 7. Required AlphaGenome notices

AlphaGenome outputs and derivatives reproduced in this paper (aggregate AVI statistics and figures) are provided under, and subject to, the AlphaGenome Output Terms of Use. By using this information, you agree to AlphaGenome Output Terms of Use found at https://deepmind.google.com/science/alphagenome/output-terms . Use of the AlphaGenome Services is subject to the AlphaGenome Services Additional Terms of Service (https://deepmind.google.com/science/alphagenome/terms, version modified 8 September 2026). AlphaGenome outputs are for non-commercial, theoretical-modelling research use only and must not be used for clinical purposes or relied on for medical or professional advice. This work is independent and is not affiliated with or endorsed by Google or Google DeepMind. The author's API credentials were personal and are not published.

## 8. Data and code availability

* **Code, pre-registered design, tests, sample lists and aggregate results:** https://github.com/Arths17/avi-clinvar-eval. A tagged release of the code is archived on Zenodo; its DOI is listed under "Related works" on this record.
* **ClinVar data** are public (NCBI); the exact file is named in Section 2.2.
* **Per-variant AlphaGenome outputs** (individual AVI scores and feature breakdowns) are **not redistributed** in the repository; they can be regenerated with the code and a personal API key.
* **GENCODE** release 46 (basic annotation, GRCh38) was used for the exploratory distance analysis (EBI GENCODE FTP, MD5 `9d4f206f82756340967a79a96db5bea6`).

## 9. Disclosures and acknowledgements

* **Independent researcher, no institutional affiliation.**
* **AI assistance:** AI tools (Claude, Anthropic) assisted with writing the code and drafting this paper. The author reviewed the work and takes responsibility for its content.
* **No clinical claims:** this work makes no clinical or diagnostic claims and is not medical advice.
* **Funding and competing interests:** No specific funding was received for this work. The author declares no competing interests (not employed or paid by Google, Google DeepMind or any related company).
* **Acknowledgements:** we thank Google DeepMind for AlphaGenome and the AlphaGenome Atlas, NCBI for ClinVar and the submitters to ClinVar, and the GENCODE consortium and EMBL-EBI for the GENCODE annotation.

## References

1. Avsec Z, *et al.* Advancing regulatory variant effect prediction with AlphaGenome. *Nature* 649, 1206–1218 (2026). doi:10.1038/s41586-025-10014-0.
2. Cheng J, *et al.* AlphaGenome Atlas: in silico mutagenesis of the entire human genome improves prioritization and interpretation of non-coding variants. *medRxiv* (preprint, not peer reviewed), doi:10.64898/2026.09.16.26363192 (2026).
3. Landrum MJ, *et al.* ClinVar: improving access to variant interpretations and supporting evidence. *Nucleic Acids Res* 46, D1062–D1067 (2018). (Cited as listed in the reference list of [1]. NCBI's ClinVar documentation asks for attribution to ClinVar as a data source. Data: https://www.ncbi.nlm.nih.gov/clinvar/, release file `clinvar_20260928.vcf.gz`.)
4. NCBI. ClinVar review status documentation, https://www.ncbi.nlm.nih.gov/clinvar/docs/review_status/ (accessed 4 October 2026).
5. AlphaGenome Services Additional Terms of Service, https://deepmind.google.com/science/alphagenome/terms (modified 8 September 2026; accessed 4 October 2026); AlphaGenome Output Terms of Use, https://deepmind.google.com/science/alphagenome/output-terms (effective 25 June 2025; accessed 4 October 2026).
6. Google DeepMind. AlphaGenome API client, https://github.com/google-deepmind/alphagenome (package version {man['alphagenome']}); documentation https://www.alphagenomedocs.com/.
7. GENCODE release 46, EMBL-EBI FTP, https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_46/ (accessed 4 October 2026).
"""

(P / "figures").mkdir(parents=True, exist_ok=True)
for f in ("roc_pooled.png", "auc_by_stratum.png", "class_counts.png"):
    shutil.copy(R / "figures" / f, P / "figures" / f)
(P / "preprint.md").write_text(TXT.strip() + "\n", encoding="utf-8")

CSS = """body{font-family:Georgia,'Times New Roman',serif;max-width:820px;margin:32px auto;padding:0 20px;line-height:1.5;font-size:11.5pt;color:#1a1a1a}
h1{font-size:19pt;line-height:1.25}h2{font-size:14pt;margin-top:1.6em;border-bottom:1px solid #ddd;padding-bottom:3px}h3{font-size:12pt}
table{border-collapse:collapse;font-size:9pt;margin:10px 0;width:100%}th,td{border:1px solid #bbb;padding:3px 6px;text-align:left}th{background:#f1f1f1}
img{max-width:80%;display:block;margin:10px auto}blockquote{border-left:4px solid #b00020;margin:12px 0;padding:4px 14px;background:#fbf3f4}
em{color:#333}code{font-size:90%}@page{margin:18mm}"""
html = markdown.markdown(TXT, extensions=["tables", "sane_lists"])
(P / "preprint.html").write_text(
    f"<!doctype html><html><head><meta charset='utf-8'><title>AVI vs ClinVar preprint</title><style>{CSS}</style></head><body>{html}</body></html>",
    encoding="utf-8")
chrome = Path("/Applications/Google Chrome.app/Contents/MacOS/Google Chrome")
if chrome.exists():
    subprocess.run([str(chrome), "--headless=new", "--disable-gpu", "--no-pdf-header-footer",
                    f"--print-to-pdf={(P / 'preprint.pdf').resolve()}", (P / "preprint.html").resolve().as_uri()],
                   check=True, capture_output=True, timeout=120)
    print("wrote paper/preprint.pdf")
print("wrote paper/preprint.md and paper/preprint.html")
import re as _re
_abs = _re.search(r"## Abstract\n\n(.*?)\n\n## 1\.", TXT, _re.S).group(1).replace("**", "")
(P / "zenodo_description.txt").write_text(
    _abs + "\n\nResearch use only. This work makes no clinical or diagnostic claims and is not medical advice. Not peer reviewed. "
    "Independent researcher, no institutional affiliation. AI tools (Claude, Anthropic) assisted with code and drafting; the author reviewed the work and takes responsibility for its content. "
    "AlphaGenome outputs reproduced here (aggregate statistics and figures) are subject to the AlphaGenome Output Terms of Use: "
    "https://deepmind.google.com/science/alphagenome/output-terms . Not affiliated with or endorsed by Google or Google DeepMind. "
    "Code and pre-registered design: https://github.com/Arths17/avi-clinvar-eval\n", encoding="utf-8")
print("wrote paper/zenodo_description.txt")
