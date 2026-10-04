import numpy as np
import pandas as pd

from avi_eval import config, sampling


def make_pool(n_genes=40, per_gene=30, seed=1):
    rng = np.random.default_rng(seed)
    rows = []
    vid = 0
    for g in range(n_genes):
        for _ in range(per_gene):
            vid += 1
            rows.append(
                {
                    "variation_id": vid, "chrom": "1", "pos": vid, "ref": "A", "alt": "G",
                    "label": int(rng.integers(0, 2)), "stars": 2,
                    "stratum": "missense", "gene": f"G{g}", "heldout_chrom": False,
                }
            )
    return pd.DataFrame(rows)


def test_deterministic():
    pool = make_pool()
    a = sampling.build_sample(pool, cohort="primary", per_class=50)
    b = sampling.build_sample(pool, cohort="primary", per_class=50)
    pd.testing.assert_frame_equal(a, b)


def test_seed_changes_sample():
    pool = make_pool()
    a = sampling.build_sample(pool, cohort="primary", per_class=50, seed=1)
    b = sampling.build_sample(pool, cohort="primary", per_class=50, seed=2)
    assert set(a["variation_id"]) != set(b["variation_id"])


def test_gene_cap_and_class_size():
    pool = make_pool()
    s = sampling.build_sample(pool, cohort="primary", per_class=50, gene_cap=3)
    assert s.groupby(["gene", "label"]).size().max() <= 3
    assert (s["label"].value_counts() <= 50).all()
    assert len(s) == 100  # plenty of eligible variants, so both classes fill up


def test_takes_all_when_cell_is_small():
    pool = make_pool(n_genes=2, per_gene=5)
    s = sampling.build_sample(pool, cohort="primary", per_class=250, gene_cap=100)
    assert len(s) == len(pool)


def test_no_gene_is_exempt_from_cap():
    pool = make_pool(n_genes=1, per_gene=40)
    pool["gene"] = None
    s = sampling.build_sample(pool, cohort="primary", per_class=1000, gene_cap=2)
    assert len(s) == 40


def test_unmapped_and_excluded_ids():
    pool = make_pool(n_genes=3, per_gene=10)
    pool.loc[:4, "stratum"] = None
    s = sampling.build_sample(pool, cohort="onestar", per_class=100, gene_cap=100, exclude_ids={10, 11})
    assert not s["variation_id"].isin(range(1, 6)).any()
    assert not s["variation_id"].isin([10, 11]).any()
    assert (s["cohort"] == "onestar").all()


def test_dry_run_is_balanced_and_reproducible():
    pool = make_pool()
    s = sampling.build_sample(pool, cohort="primary", per_class=50)
    d1, d2 = sampling.pick_dry_run(s), sampling.pick_dry_run(s)
    pd.testing.assert_frame_equal(d1, d2)
    assert len(d1) == config.DRY_RUN_N
    assert d1["label"].sum() == config.DRY_RUN_N // 2
