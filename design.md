# Study design (pre-registration) — DRAFT v0.1, awaiting author approval

**Working title:** How well does the AlphaGenome Variant Impact (AVI) score separate pathogenic from benign human genetic variants? An independent evaluation on ClinVar.
**Author:** Atharv Ranjan, independent researcher.
**Status:** No AVI scores have been retrieved or analysed. This file is frozen once approved; any later change goes in a dated "Deviations" section at the bottom, never as a silent edit.

> Research use only. Nothing here is a clinical or diagnostic claim or medical advice.

---

## 0. Facts this design rests on (all read from primary sources on 2026-10-04)

| Fact | Source |
|---|---|
| Atlas covers all ~9 billion SNVs on **GRCh38/hg38**; Atlas website has no precomputed indel AVI | Atlas FAQ |
| Python access: `atlas.create(api_key)` → `AtlasClient.query_variants(...)`; scorers `AVI_SCORE` and `AVI_SCORE_FEATURE_IMPORTANCE` | `alphagenome` repo (`atlas.py`); DeepMind `science-skills` SKILL.md |
| No numeric quota is published; "thousands" of on-demand predictions is the intended scale, >1 million unsuitable | alphagenomedocs.com, repo README |
| AVI = 18-feature supervised model, trained on gnomAD-style proxy labels (filtering AF above vs below 0.1%), **not on ClinVar labels** | Atlas preprint (Cheng et al. 2026, medRxiv) |
| AVI training variants come from chr 1,4,7,8,10,13,15; validation from chr 2,5,11,14,17,20,22,X | Atlas preprint, Methods |
| Atlas authors benchmarked AVI on ClinVar (June 15, 2025 release, SNVs and indels, **no star filter**) | Atlas preprint, Methods |
| AVI inputs include AlphaMissense and loss-of-function flags derived partly from "ClinVar nonsense / initiator_codon_variant" descriptors | Atlas preprint, Methods |
| ClinVar: monthly archived releases (first Thursday), GRCh38 VCF available; attribution requested; not for diagnostic use | NCBI ClinVar docs |

---

## 1. Variant source

- **ClinVar release:** the monthly GRCh38 VCF released **2026-10-01** (first Thursday of October 2026). If that file is unavailable, use the most recent archived monthly release and record its date. The exact file name and its SHA-256 are recorded in `results/run_manifest.json`.
- **Variant type:** single-nucleotide variants only (one reference base, one alternate base). Reason: Atlas AVI is precomputed for SNVs; indels are out of scope.
- **Genome build:** GRCh38 (matches Atlas). No liftover is used.
- **Chromosomes:** autosomes and X. Y and mitochondrial variants are excluded (the Atlas authors also excluded them).
- **One row per ClinVar VariationID.**

## 2. Labels

*ClinVar review status:* ClinVar gives each record a 0–4 star rating for how well-supported the classification is (1 star = a single submitter with criteria; 2 stars = several submitters who agree; 3 = expert panel; 4 = practice guideline).

- **Positive (pathogenic):** aggregate classification exactly `Pathogenic`, `Likely_pathogenic`, or `Pathogenic/Likely_pathogenic`.
- **Negative (benign):** exactly `Benign`, `Likely_benign`, or `Benign/Likely_benign`.
- **Excluded:** uncertain significance, `Conflicting_classifications_of_pathogenicity`, risk factor, drug response, low-penetrance terms, and any record whose classification string also contains a term outside the lists above.
- **Primary review-status threshold: ≥2 stars** (`criteria_provided,_multiple_submitters,_no_conflicts`, `reviewed_by_expert_panel`, `practice_guideline`).
  - *Why:* single-submitter labels are the noisiest. Requiring agreement or expert review gives cleaner labels.
  - *Cost:* the pool shifts toward well-studied genes, and benign labels from large submitters dominate. Both are listed in Limitations.
- **Sensitivity analysis:** repeat the primary analysis with ≥1 star (the Atlas authors used no filter, so this is the most comparable setting to theirs).

## 3. Consequence categories (strata)

*Molecular consequence:* a standard label (Sequence Ontology term) for what a variant does to a gene, such as "missense" (changes one amino acid) or "synonymous" (changes the DNA but not the amino acid). We take it from ClinVar's own `MC` field, so no extra download is needed.

Planned strata (mapped from the SO terms in `MC`; if a variant has several terms across transcripts, the first matching row of this priority list wins):

1. nonsense / stop-gained, start-lost, stop-lost
2. canonical splice site (splice donor / splice acceptor)
3. missense
4. synonymous
5. splice region (if separable in `MC`; otherwise merged into intronic)
6. intronic
7. 5′ UTR
8. 3′ UTR
9. other non-coding (ncRNA, up/downstream, etc.)

**Gate:** after the ClinVar VCF is parsed, I inspect only the counts of `MC` terms and finalise this mapping *before any API call*. The final mapping is committed with its date. Labels/scores are never used to adjust it.

Roll-up groups used in reporting: **coding** (1, 3, 4), **splice** (2, 5), **non-coding** (6–9).

## 4. Sampling

- **Cap per gene:** at most **10 variants per gene per class per stratum**, chosen at random. Reason: BRCA1/BRCA2/MLH1-type genes carry thousands of ClinVar variants and would dominate otherwise.
- **Per stratum:** up to **250 pathogenic and 250 benign** variants, sampled uniformly at random from what remains. If a class has fewer, take all of them.
- **Strata with fewer than 30 variants in either class** are described (counts, scores) but get **no AUC**.
- **Total budget:** at most ~4,000 variants (≤ 9 strata × 500). **Hard cap of 6,000 API queries** including retries and dry run.
- **Why this size:** the Atlas API is intended for thousands of queries, and 4,000 is far below the 1-million ceiling. A rough Hanley–McNeil calculation for AUC≈0.85 with 250 + 250 variants gives a 95% CI half-width of about ±0.03 (approximate, to be re-checked in code). That is enough to compare categories to within a few hundredths, not to detect tiny differences.
- **Seed:** `20261004` (fixed for numpy, Python `random`, and bootstrap).
- **Attrition log:** every filter step records how many variants remain (`results/attrition.csv`), in the order: raw VCF rows → SNVs → GRCh38 chr1–22,X → label sets → review ≥2 stars → mapped stratum → gene cap → final sample.
- **Pre-sampled held-out-chromosome flag:** each variant is tagged `heldout_chrom = chr ∈ {3,6,9,12,16,18,19,21}` (the autosomes not used for AVI training or validation, per the Atlas preprint). Not used for sampling, only for the sensitivity analysis in §6.

## 5. Scores and primary metrics

- **Score:** raw `AVI_SCORE` (higher = predicted more impactful). Using the raw value avoids ties from the capped PHRED scale. PHRED and quantile are stored too (*PHRED:* a log scale where 10 = top 10% of all genome variants, 20 = top 1%, and so on).
- **ROC-AUC** (*the chance that a randomly chosen pathogenic variant gets a higher score than a randomly chosen benign one; 0.5 = coin flip, 1.0 = perfect*) and **PR-AUC** (*average precision; how well the top-ranked variants are actually pathogenic; its coin-flip level equals the share of pathogenic variants, so it must always be read next to that prevalence*).
- **Confidence intervals:** 95% percentile bootstrap, 2,000 replicates, **resampling genes (not variants)** because variants in one gene are not independent. Variants without a gene symbol form their own singleton clusters.
- **Primary result:** pooled ROC-AUC and PR-AUC over the whole sample, **plus** the macro-average of per-stratum ROC-AUCs (each eligible stratum counts equally). The pooled number is inflated by consequence mix (see §7), so conclusions are drawn from the stratum-level numbers and from the comparison with the baseline.
- **Secondary:** per-stratum and per-roll-up ROC-AUC/PR-AUC with CIs; sensitivity (≥1 star); held-out-chromosome subset; fixed-threshold sensitivity/specificity at PHRED ≥ 10, 20, 30 (thresholds fixed now, never tuned).
- **Direction is fixed:** higher AVI = pathogenic. AUCs below 0.5 are reported as they are, not flipped.
- **No model is trained or fitted on any AlphaGenome output** (including recalibration or logistic regression on AVI). Only rank-based metrics are computed. This keeps the project clear of the ML-training restriction.

## 6. Baselines and ablations

1. **Consequence-only baseline (B1):** a fixed ordinal "severity" score by stratum, set now and not tuned: nonsense/start/stop-lost = 5 > canonical splice = 4 > missense = 3 > splice region = 2 > synonymous = 1 > intronic = UTR = other non-coding = 0 (ties allowed). It uses no data fitting. It answers "does AVI add anything beyond knowing the consequence type?" and is meaningful only for the pooled analysis (within a stratum it is constant, AUC 0.5 by construction).
   - *Caveat:* the consequence label comes from ClinVar annotations, and AVI's own loss-of-function inputs also draw on ClinVar descriptors, so B1 and AVI share some information. This is stated in the paper.
2. **AVI decomposition (A1):** `AVI_SCORE_FEATURE_IMPORTANCE` gives each variant's 18 SHAP contributions (*SHAP: a method that splits one prediction into additive contributions from each input feature*). I compute (i) the AlphaMissense contribution alone and (ii) the sum of the ten AlphaGenome-derived contributions alone, and report their AUCs per stratum. Purpose: show how much of AVI's performance comes from the AlphaMissense input (relevant to circularity) versus the sequence model. The feature breakdown is not part of the "AVI Score" in the terms, so it is treated as ordinary restricted Output (non-commercial; notice attached).
3. **Held-out-chromosome analysis:** repeat the primary analysis on `heldout_chrom` variants only. This removes AVI's own training/validation chromosomes. It does **not** remove possible chromosome overlap with the underlying AlphaGenome or AlphaMissense training, which I cannot verify.
4. **No external baseline (CADD, REVEL, etc.)** is used, because it would need large downloads. So the paper will **not** claim AVI is better or worse than other methods, only how it performs here relative to B1 and to its own components.

## 7. Confounds and checks (all reported, none optional)

- **Class imbalance across consequence types:** a table of positive/negative counts per stratum. Nonsense is almost all pathogenic and synonymous almost all benign, so a pooled AUC mostly reflects consequence mix. Hence the stratum-level focus and B1.
- **Gene clustering:** gene-level bootstrap (above); plus a leave-top-5-genes-out recomputation.
- **ClinVar circularity (the main risk).** Known from sources: (a) AVI training labels are not ClinVar, but its loss-of-function input flags use ClinVar descriptors; (b) the Atlas authors themselves evaluate on ClinVar, so performance there is the headline claim being replicated, not an independent test; (c) AlphaMissense is an input, and I have not verified whether ClinVar was used to train or calibrate it; (d) the AlphaGenome Nature paper lists ClinVar among training/evaluation data sources in a supplement I have not read. **If leakage cannot be excluded, the abstract and limitations must say so prominently** (it cannot be excluded here).
- **Temporal check (optional, only if the author approves a second ClinVar download):** compare against the June 2025 archived release and report performance on variants absent from it (i.e., classified after the Atlas authors' evaluation snapshot).
- **Label noise / ClinVar bias:** ClinVar over-represents well-studied genes and European-ancestry cohorts; pathogenic labels from single labs may be wrong. Discussed, not corrected.
- **Sanity gate (Stage 3):** if pooled ROC-AUC > 0.97 or < 0.55, or any stratum AUC is implausibly extreme, stop and look for bugs (wrong allele, coordinate off-by-one, label flip, duplicate rows, leakage) *before* interpreting.

## 8. API usage, caching, secrets

- Key read from env var `ALPHAGENOME_API_KEY` (also accepts `ALPHA_GENOME_API_KEY`) or Kaggle Secrets. Never logged, printed, or committed.
- Every response cached to disk keyed by `chrom:pos:ref>alt`; reruns make zero new calls.
- Conservative concurrency (4 workers), exponential backoff with jitter on timeouts/transient errors, immediate stop on permission errors.
- Dry run: 20 variants.

## 9. What gets published

- **Public:** code, this design, the list of sampled ClinVar VariationIDs with coordinates and labels (ClinVar data are unrestricted), aggregate metrics, bootstrap summaries, plots.
- **Not public (until the author has re-checked the terms in Stage 5 and decided):** per-variant AVI values and feature breakdowns. They stay in a git-ignored local cache.
- **Notices in the repo and paper:** the Output Terms notice ("By using this information, you agree to AlphaGenome Output Terms of Use found at https://deepmind.google.com/science/alphagenome/output-terms"), citations to Avsec et al. 2026 (Nature) and Cheng et al. 2026 (medRxiv), ClinVar attribution (Landrum et al.), non-commercial and no-clinical-use statements.

## 10. Known limitations to carry into the paper

Circularity risk (§7); single ClinVar release and single genome build; SNVs only; ≥2-star filter biases toward well-studied genes and common variant types; per-gene results not reported because of the per-gene cap; no external baseline; small sample per stratum, so wide CIs for rare categories; independent evaluation by one researcher, with AI assistance (Claude) for code and drafting.

## Deviations log

Design approved by the author on 2026-10-04. Entries below were made after approval and before any AVI score was retrieved.

1. **2026-10-04 — ClinVar release file.** §1 assumed a monthly VCF released 2026-10-01. NCBI's GRCh38 VCF directory shows dated files at roughly weekly cadence, and the newest was `clinvar_20260928.vcf.gz` (~189 MB); no 2026-10-01 file existed. Per the pre-registered fallback in §1, the primary release is **`clinvar_20260928.vcf.gz`**. Its MD5 and SHA-256 are written to `results/run_manifest.json`.
2. **2026-10-04 — Scope of the ≥1-star sensitivity analysis (§2).** Repeating the full primary sampling at ≥1 star could need up to ~4,500 extra API calls, which would break the 6,000-call cap in §4. Instead the ≥1-star analysis uses the primary sample **plus a small supplementary sample of exactly-1-star variants** (up to 60 per class per stratum, same gene cap and seed logic). The "≥1 star" set is the union of the two. Its class and stratum mix therefore differs from the primary sample, and the paper will say so.
4. **2026-10-04 — Stratum-mapping gate (§3), decided from term counts only, before any API call.** On the real `clinvar_20260928` file: (a) ClinVar's `MC` field contains **no `splice_region_variant` term**, so the `splice_region` stratum is **empty**; per §3 ("otherwise merged into intronic") such variants sit in `intronic` (or `missense`/`synonymous` when they carry those terms). The "splice" roll-up is therefore canonical splice donor/acceptor only. (b) **2,067** otherwise-eligible variants (≥2 stars) have **no `MC` value at all**, and 85 benign variants carry only `no_sequence_alteration`; both groups are unmapped and excluded (logged in `results/attrition.csv`). No other term was unmapped.
5. **2026-10-04 — Expected class imbalance in non-coding strata.** After the gene cap, pathogenic variants available at ≥2 stars were: `other_noncoding` 17, `utr3` 22, `utr5` 48 (benign: 250 each). Under the pre-registered rule (§4: ≥30 per class for an AUC), `other_noncoding` and `utr3` will be **described but get no AUC**; `utr5` is eligible but will have a wide CI. The "non-coding" roll-up (which pools them with intronic) is still reported. Sample sizes are not increased, to keep within the API budget. Real sample: 3,228 primary + 887 one-star = 4,115 variants.
6. **2026-10-04 — Sanity gate triggered; post-hoc diagnostics (EXPLORATORY, not pre-registered).** After the full run, the §7 gate flagged ROC-AUC > 0.97 in `synonymous` (0.975) and `intronic` (0.980). The primary analysis and all its parameters were NOT changed. Diagnostics (script `scripts/diagnostics.py`, output `results/diagnostics*.{json,csv}`) were run to look for bugs and leakage: (a) integrity — unique IDs/keys, independent scikit-learn AUCs reproduce the pipeline's values, and 300 randomly chosen sampled variants re-read from the raw VCF show 0 coordinate/label mismatches; (b) a within-stratum label-shuffle null gives macro-AUC ≈ 0.50; (c) the signal in these two strata is carried by the splicing feature (single-feature AUC 0.989 synonymous, 0.976 intronic), while the AlphaMissense contribution is ≈ 0.5 there; (d) population-frequency shortcut: ClinVar benign variants far more often carry an ExAC/1000G/ESP allele frequency (62–96%) than pathogenic ones (8–29%), but AUC among variants with *no* such frequency stays high (intronic 0.979, synonymous 0.987); (e) restricted to chromosomes AVI did not train or validate on, AUCs are essentially unchanged (intronic 0.980, synonymous 0.977). Conclusion recorded: no bug was found and no leakage via AVI's own training chromosomes or a frequency shortcut was found; the high values are reported as observed, with the caveat that ClinVar pathogenic synonymous/intronic variants are likely enriched for splice-disrupting changes near exons (not tested: no exon annotation was downloaded), that independence from ClinVar cannot be guaranteed (§7), and that the small no-frequency benign groups (47 intronic, 96 synonymous) make those sub-estimates imprecise.
7. **2026-10-04 — Splice-proximity diagnostic (EXPLORATORY, not pre-registered).** With the author's approval, `gencode.v46.basic.annotation.gtf.gz` (GENCODE release 46, GRCh38, ~29 MB, MD5 `9d4f206f82756340967a79a96db5bea6` verified against the EBI `MD5SUMS` file) was downloaded to test the untested explanation in deviation 6. `scripts/splice_distance.py` computes, for each sampled variant, the distance to the nearest exon edge in that annotation (bands fixed before looking: 0–2, 3–10, 11–50, 51–200, >200 bp) and reports AUC where each class has ≥10 variants (looser than the pre-registered 30). Outputs: `results/diagnostics_splice_distance*.csv`. Findings: (a) **composition** — pathogenic synonymous variants sit overwhelmingly at exon edges (167/250 within 0–2 bp vs 8/250 benign) and pathogenic intronic variants at 3–10 bp (163/250 vs 67/250 benign), i.e. in these strata the label is strongly associated with proximity to exon boundaries, as hypothesised; (b) **within-band discrimination persists** — AVI AUC intronic 0.996 (3–10 bp), 0.961 (11–50), 0.937 (51–200), 0.961 (>200); synonymous 0.973 (11–50), 0.920 (51–200); gene-clustered bootstrap 95% CIs (2,000 replicates, same method as the main analysis) were then added: intronic 11–50 bp 0.961 [0.924, 0.989], 51–200 bp 0.937 [0.853, 0.993], >200 bp 0.961 [0.901, 1.000]; synonymous 11–50 bp 0.973 [0.927, 0.998], 51–200 bp 0.920 [0.819, 0.989]; several bands remain small (e.g. 11 pathogenic intronic variants at 51–200 bp), so these remain indicative. Sanity check: no intronic variant has distance 0–2 because distances 1–2 are canonical splice sites, which are in a different stratum. Caveat: nearest edge is taken over all exons in the GENCODE basic set, regardless of gene or transcript.
3. **2026-10-04 — Extra definition detail.** The gene used for the cap and the bootstrap is the **first gene symbol** listed in the ClinVar `GENEINFO` field. Variants with no gene symbol are each treated as their own cluster and are exempt from the cap.
