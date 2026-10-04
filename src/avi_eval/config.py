"""Fixed study parameters. Every value here is pre-registered in design.md.

Do not change a value after scores have been retrieved; record any change in
the "Deviations log" of design.md instead.
"""
from __future__ import annotations

from pathlib import Path

SEED = 20261004

# --- Paths (relative to the repository root, i.e. the current directory) -----
ROOT = Path(".")
DATA_DIR = ROOT / "data"
RAW_DIR = DATA_DIR / "raw"
RESULTS_DIR = ROOT / "results"
CACHE_DIR = ROOT / "cache"

# --- ClinVar source (design.md section 1, deviation 1) -----------------------
CLINVAR_RELEASE = "20260928"
CLINVAR_BASE_URL = "https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38"

# --- Variants -----------------------------------------------------------------
CHROMS = tuple(str(i) for i in range(1, 23)) + ("X",)
# Autosomes used for neither AVI training nor validation (Atlas preprint, Methods).
HELDOUT_CHROMS = frozenset({"3", "6", "9", "12", "16", "18", "19", "21"})

# --- Labels (design.md section 2) ---------------------------------------------
POSITIVE_TERMS = frozenset({"pathogenic", "likely_pathogenic", "pathogenic/likely_pathogenic"})
NEGATIVE_TERMS = frozenset({"benign", "likely_benign", "benign/likely_benign"})

# ClinVar review-status string -> star rating (NCBI "review status" documentation).
# Keys have commas removed, because the VCF writes e.g.
# "criteria_provided,_multiple_submitters,_no_conflicts".
REVIEW_STARS = {
    "practice_guideline": 4,
    "reviewed_by_expert_panel": 3,
    "criteria_provided_multiple_submitters_no_conflicts": 2,
    "criteria_provided_conflicting_classifications": 1,
    "criteria_provided_conflicting_interpretations": 1,  # older wording
    "criteria_provided_single_submitter": 1,
    "no_assertion_criteria_provided": 0,
    "no_classification_provided": 0,
    "no_classification_for_the_single_variant": 0,
    "no_classifications_from_unflagged_records": 0,
    "no_assertion_provided": 0,
}
PRIMARY_MIN_STARS = 2

# --- Strata (design.md section 3) ---------------------------------------------
# Priority order: the first stratum whose terms intersect the variant's MC terms wins.
# Term lists are lower-case Sequence Ontology names as written in the ClinVar MC field.
# GATE: after `prepare`, inspect data/mc_term_counts.csv and finalise this mapping
# BEFORE any API call. Never adjust it using scores.
STRATUM_TERMS: dict[str, frozenset[str]] = {
    "nonsense_start_stop": frozenset({"nonsense", "stop_gained", "stop_lost", "initiator_codon_variant", "start_lost"}),
    "splice_canonical": frozenset({"splice_donor_variant", "splice_acceptor_variant"}),
    "missense": frozenset({"missense_variant"}),
    "synonymous": frozenset({"synonymous_variant"}),
    "splice_region": frozenset({"splice_region_variant"}),
    "intronic": frozenset({"intron_variant"}),
    "utr5": frozenset({"5_prime_utr_variant"}),
    "utr3": frozenset({"3_prime_utr_variant"}),
    "other_noncoding": frozenset(
        {
            "non-coding_transcript_variant",
            "non_coding_transcript_variant",
            "non-coding_transcript_exon_variant",
            "non_coding_transcript_exon_variant",
            "upstream_transcript_variant",
            "downstream_transcript_variant",
            "genic_upstream_transcript_variant",
            "genic_downstream_transcript_variant",
            "2kb_upstream_variant",
            "500b_downstream_variant",
            "upstream_gene_variant",
            "downstream_gene_variant",
        }
    ),
}
STRATA = tuple(STRATUM_TERMS)

ROLLUPS = {
    "coding": ("nonsense_start_stop", "missense", "synonymous"),
    "splice": ("splice_canonical", "splice_region"),
    "non_coding": ("intronic", "utr5", "utr3", "other_noncoding"),
}

# Baseline B1: fixed ordinal severity (design.md section 6). Never tuned.
SEVERITY = {
    "nonsense_start_stop": 5,
    "splice_canonical": 4,
    "missense": 3,
    "splice_region": 2,
    "synonymous": 1,
    "intronic": 0,
    "utr5": 0,
    "utr3": 0,
    "other_noncoding": 0,
}

# --- Sampling (design.md section 4) -------------------------------------------
PER_CLASS_PER_STRATUM = 250
GENE_CAP = 10  # per gene, per class, per stratum
ONESTAR_PER_CLASS_PER_STRATUM = 60  # supplementary sample, deviation 2
MIN_PER_CLASS_FOR_AUC = 30
MAX_SAMPLE_TOTAL = 4500  # 9 strata x 500
DRY_RUN_N = 20

# --- Metrics (design.md section 5) --------------------------------------------
N_BOOT = 2000
PHRED_THRESHOLDS = (10, 20, 30)  # fixed; never tuned
SANITY_AUC_HIGH = 0.97
SANITY_AUC_LOW = 0.55

# --- API usage (design.md section 8) ------------------------------------------
MAX_API_CALLS = 6000  # hard cap, including retries and dry run
WORKERS = 4
MAX_ATTEMPTS = 5
BACKOFF_BASE_SECONDS = 2.0
API_KEY_ENV_NAMES = ("ALPHAGENOME_API_KEY", "ALPHA_GENOME_API_KEY")
REQUESTED_SCORERS = ("AVI_SCORE", "AVI_SCORE_FEATURE_IMPORTANCE")

# AVI feature names (DeepMind science-skills `AviFeature` enum). The first ten are
# derived from AlphaGenome predictions; see design.md section 6 (A1).
ALPHAGENOME_FEATURES = (
    "MERGED_SPLICING",
    "MAX_ABS_RNA_SEQ",
    "MAX_ABS_ATAC",
    "MAX_ABS_DNASE",
    "MAX_ABS_CHIP_TF",
    "MAX_ABS_CHIP_HISTONE",
    "MAX_ABS_CAGE",
    "MAX_ABS_PROCAP",
    "MAX_ABS_POLYADENYLATION",
    "MAX_ABS_CONTACT_MAPS",
)
ALPHAMISSENSE_FEATURE = "ALPHAMISSENSE"
OTHER_FEATURES = (
    "CACTUS_241_WAY",
    "PROTEIN_TERMINATION",
    "START_LOST",
    "STOP_LOST",
    "PHASTCONS_470_WAY",
    "IS_INSERTION",
    "IS_DELETION",
)
