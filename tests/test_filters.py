import pytest

from avi_eval import clinvar, config


@pytest.mark.parametrize(
    "sig,expected",
    [
        ("Pathogenic", 1),
        ("Likely_pathogenic", 1),
        ("Pathogenic/Likely_pathogenic", 1),
        ("Benign", 0),
        ("Likely_benign", 0),
        ("Benign/Likely_benign", 0),
        ("Uncertain_significance", None),
        ("Conflicting_classifications_of_pathogenicity", None),
        ("Pathogenic|risk_factor", None),
        ("Pathogenic|drug_response", None),
        ("Pathogenic,_low_penetrance", None),
        ("Pathogenic|Benign", None),
        ("", None),
        (None, None),
    ],
)
def test_label_from_clnsig(sig, expected):
    assert clinvar.label_from_clnsig(sig) == expected


@pytest.mark.parametrize(
    "text,stars",
    [
        ("practice_guideline", 4),
        ("reviewed_by_expert_panel", 3),
        ("criteria_provided,_multiple_submitters,_no_conflicts", 2),
        ("criteria_provided_multiple_submitters_no_conflicts", 2),
        ("criteria_provided,_single_submitter", 1),
        ("criteria_provided,_conflicting_classifications", 1),
        ("no_assertion_criteria_provided", 0),
        ("something_new", -1),
        (None, -1),
    ],
)
def test_stars(text, stars):
    assert clinvar.stars_from_revstat(text) == stars


def test_is_snv():
    assert clinvar.is_snv("A", "G")
    assert not clinvar.is_snv("A", "A")
    assert not clinvar.is_snv("A", "GA")
    assert not clinvar.is_snv("AT", "A")
    assert not clinvar.is_snv("N", "A")


def test_normalise_chrom():
    assert clinvar.normalise_chrom("chr7") == "7"
    assert clinvar.normalise_chrom("X") == "X"


def test_parse_mc_and_first_gene():
    assert clinvar.parse_mc("SO:0001583|missense_variant,SO:0001627|intron_variant") == [
        "missense_variant",
        "intron_variant",
    ]
    assert clinvar.parse_mc(None) == []
    assert clinvar.first_gene("BRCA1:672|BRCA1-AS:1234") == "BRCA1"
    assert clinvar.first_gene(None) is None
    assert clinvar.first_gene(":123") is None


def test_stratum_priority():
    # missense (3) beats splice_region (5) and intron (6)
    assert clinvar.assign_stratum(["splice_region_variant", "missense_variant"]) == "missense"
    assert clinvar.assign_stratum(["intron_variant", "splice_region_variant"]) == "splice_region"
    assert clinvar.assign_stratum(["splice_donor_variant", "intron_variant"]) == "splice_canonical"
    assert clinvar.assign_stratum(["nonsense"]) == "nonsense_start_stop"
    assert clinvar.assign_stratum(["coding_sequence_variant"]) is None
    assert clinvar.assign_stratum([]) is None


def test_every_stratum_has_severity_and_rollup():
    assert set(config.SEVERITY) == set(config.STRATA)
    rolled = [s for members in config.ROLLUPS.values() for s in members]
    assert sorted(rolled) == sorted(config.STRATA)


def test_check_header():
    ok = [f"##INFO=<ID={i},Number=.,Type=String>" for i in clinvar.REQUIRED_INFO_IDS]
    clinvar.check_header(ok)
    with pytest.raises(ValueError, match="MC"):
        clinvar.check_header([x for x in ok if "ID=MC," not in x])


def test_load_clinvar_counts_and_attrition(tiny_vcf):
    att = clinvar.Attrition()
    df, unknown = clinvar.load_clinvar(tiny_vcf, att)
    steps = {r["step"]: r["n_remaining"] for r in att.rows}
    assert steps["raw VCF data rows"] == 10
    assert steps["single-nucleotide variants"] == 9      # one indel removed
    assert steps["on chr1-22 or X"] == 8                 # MT removed
    assert steps["clean pathogenic-type or benign-type label"] == 6  # ids 1,2,3,7,9,10
    assert sorted(df["variation_id"]) == [1, 2, 3, 9]
    assert unknown["some_new_status"] == 1
    assert set(df["label"]) == {0, 1}
    assert df.set_index("variation_id").loc[9, "stars"] == 3
    assert df.set_index("variation_id").loc[3, "stars"] == 1


def test_add_strata_and_heldout(tiny_vcf):
    df, _ = clinvar.load_clinvar(tiny_vcf, clinvar.Attrition())
    df = clinvar.add_strata(df).set_index("variation_id")
    assert df.loc[1, "stratum"] == "missense"
    assert df.loc[2, "stratum"] == "synonymous"
    assert bool(df.loc[1, "heldout_chrom"]) is False  # chr1 was an AVI training chromosome
    assert config.HELDOUT_CHROMS.isdisjoint({"1", "4", "7", "8", "10", "13", "15", "2", "5", "11", "14", "17", "20", "22", "X"})
