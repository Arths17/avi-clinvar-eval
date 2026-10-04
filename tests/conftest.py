import gzip
import textwrap

import pytest

HEADER = textwrap.dedent(
    """\
    ##fileformat=VCFv4.1
    ##INFO=<ID=CLNSIG,Number=.,Type=String,Description="x">
    ##INFO=<ID=CLNREVSTAT,Number=.,Type=String,Description="x">
    ##INFO=<ID=MC,Number=.,Type=String,Description="x">
    ##INFO=<ID=GENEINFO,Number=1,Type=String,Description="x">
    ##INFO=<ID=CLNVC,Number=1,Type=String,Description="x">
    #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO
    """
)


def vcf_line(chrom, pos, vid, ref, alt, sig, rev, mc="SO:0001583|missense_variant", gene="GENE1:1"):
    info = f"CLNSIG={sig};CLNREVSTAT={rev};MC={mc};GENEINFO={gene};CLNVC=single_nucleotide_variant"
    return f"{chrom}\t{pos}\t{vid}\t{ref}\t{alt}\t.\t.\t{info}\n"


@pytest.fixture
def tiny_vcf(tmp_path):
    two = "criteria_provided,_multiple_submitters,_no_conflicts"
    one = "criteria_provided,_single_submitter"
    zero = "no_assertion_criteria_provided"
    lines = [
        vcf_line("1", 100, 1, "A", "G", "Pathogenic", two),                       # keep (primary)
        vcf_line("1", 200, 2, "C", "T", "Benign", two, mc="SO:0001819|synonymous_variant"),  # keep
        vcf_line("2", 300, 3, "G", "A", "Likely_pathogenic", one),               # keep, 1 star
        vcf_line("2", 400, 4, "G", "A", "Uncertain_significance", two),          # drop: VUS
        vcf_line("2", 500, 5, "G", "GA", "Pathogenic", two),                     # drop: indel
        vcf_line("MT", 600, 6, "A", "G", "Pathogenic", two),                     # drop: chrom
        vcf_line("3", 700, 7, "A", "C", "Pathogenic", zero),                     # drop: 0 stars
        vcf_line("3", 800, 8, "A", "C", "Conflicting_classifications_of_pathogenicity", two),  # drop
        vcf_line("X", 900, 9, "T", "C", "Benign/Likely_benign", "reviewed_by_expert_panel"),  # keep
        vcf_line("4", 950, 10, "T", "C", "Pathogenic", "some_new_status"),       # drop: unknown status
    ]
    path = tmp_path / "tiny.vcf.gz"
    with gzip.open(path, "wt") as fh:
        fh.write(HEADER)
        fh.writelines(lines)
    return path
