import numpy as np
import pandas as pd
import pytest
from sklearn.metrics import average_precision_score, roc_auc_score

from avi_eval import analysis, config, metrics, scoring


def test_roc_auc_matches_sklearn_including_ties():
    rng = np.random.default_rng(0)
    y = rng.integers(0, 2, 300)
    s = rng.integers(0, 10, 300).astype(float)  # many ties
    assert metrics.roc_auc(y, s) == pytest.approx(roc_auc_score(y, s))


def test_roc_auc_known_values():
    assert metrics.roc_auc([0, 0, 1, 1], [0.1, 0.2, 0.8, 0.9]) == 1.0
    assert metrics.roc_auc([0, 0, 1, 1], [0.9, 0.8, 0.2, 0.1]) == 0.0
    assert metrics.roc_auc([0, 1], [0.5, 0.5]) == 0.5
    assert np.isnan(metrics.roc_auc([1, 1, 1], [1, 2, 3]))


def test_average_precision():
    y = np.array([0, 1, 0, 1, 1, 0])
    s = np.array([0.1, 0.9, 0.3, 0.7, 0.6, 0.2])
    assert metrics.avg_precision(y, s) == pytest.approx(average_precision_score(y, s))
    assert np.isnan(metrics.avg_precision([0, 0], [1, 2]))


def test_severity_baseline_is_fixed_and_ordinal():
    b = metrics.severity_baseline(["nonsense_start_stop", "missense", "synonymous", "intronic"])
    assert list(b) == [5, 3, 1, 0]


def synthetic(n_genes=60, per_gene=12, signal=1.5, seed=3):
    rng = np.random.default_rng(seed)
    rows, vid = [], 0
    for g in range(n_genes):
        for _ in range(per_gene):
            vid += 1
            label = int(rng.integers(0, 2))
            rows.append(
                {
                    "variation_id": vid, "gene": f"G{g}", "label": label,
                    "stratum": "missense" if vid % 2 else "synonymous",
                    "score": rng.normal(signal * label, 1.0),
                }
            )
    return pd.DataFrame(rows)


def test_bootstrap_reproducible_and_brackets_point_estimate():
    df = synthetic()
    groups = {"pooled": np.ones(len(df), bool)}
    a = metrics.bootstrap(df, ["score"], groups, n_boot=200, seed=5)
    b = metrics.bootstrap(df, ["score"], groups, n_boot=200, seed=5)
    np.testing.assert_array_equal(a["boot_roc"][("score", "pooled")], b["boot_roc"][("score", "pooled")])
    point = a["point"][("score", "pooled")][0]
    lo, hi = metrics.percentile_ci(a["boot_roc"][("score", "pooled")])
    assert lo < point < hi
    assert 0.7 < point < 0.95  # sanity of the synthetic signal, not a study result


def test_cluster_bootstrap_wider_than_row_bootstrap_when_clustered():
    # Make labels and scores strongly gene-dependent so clustering matters.
    rng = np.random.default_rng(9)
    rows, vid = [], 0
    for g in range(30):
        gene_label = g % 2
        gene_shift = rng.normal(0, 1)
        for _ in range(20):
            vid += 1
            rows.append({"variation_id": vid, "gene": f"G{g}", "label": gene_label,
                         "stratum": "missense", "score": gene_shift + rng.normal(0, 0.3)})
    df = pd.DataFrame(rows)
    groups = {"pooled": np.ones(len(df), bool)}
    clustered = metrics.bootstrap(df, ["score"], groups, n_boot=300, seed=1)
    df2 = df.assign(gene=None)  # every row its own cluster -> ordinary bootstrap
    plain = metrics.bootstrap(df2, ["score"], groups, n_boot=300, seed=1)
    w = lambda r: np.subtract(*metrics.percentile_ci(r["boot_roc"][("score", "pooled")])[::-1])
    assert w(clustered) > w(plain)


def test_eligibility_threshold():
    y = np.array([1] * 29 + [0] * 40)
    assert not metrics.eligible(y)
    assert metrics.eligible(np.array([1] * 30 + [0] * 30))


def test_leave_top_genes_out_removes_most_frequent():
    df = synthetic(n_genes=10, per_gene=10)
    df.loc[:49, "gene"] = "BIG"
    out = metrics.leave_top_genes_out(df, "score", k=1)
    assert out["removed_genes"] == ["BIG"]
    assert out["n_remaining"] == 50


def test_phred_conversion():
    assert scoring.phred_from_cdf_quantile(0.9) == pytest.approx(10.0)
    assert scoring.phred_from_cdf_quantile(0.99) == pytest.approx(20.0)
    assert np.isfinite(scoring.phred_from_cdf_quantile(1.0))  # clipped, no inf


# ---- parsing API responses and the feature decomposition -------------------- #
def fake_adata(values, names=None, quantile=None):
    import anndata

    X = np.array(values, dtype=np.float32).reshape(1, -1)
    var = pd.DataFrame(index=[str(i) for i in range(X.shape[1])])
    if names is not None:
        var["name"] = names
    layers = None if quantile is None else {"quantiles": np.array([[quantile]], dtype=np.float32)}
    return anndata.AnnData(X=X, var=var, layers=layers)


def all_feature_names():
    return list(config.ALPHAGENOME_FEATURES) + [config.ALPHAMISSENSE_FEATURE] + list(config.OTHER_FEATURES)


def test_parse_avi_response_ok():
    names = all_feature_names()
    resp = {
        "AVI_SCORE": fake_adata([3.5], quantile=0.99),
        "AVI_SCORE_FEATURE_IMPORTANCE": fake_adata(np.arange(len(names)), names),
    }
    rec = scoring.parse_avi_response(resp)
    assert rec["avi_raw"] == pytest.approx(3.5)
    assert rec["avi_phred"] == pytest.approx(20.0, abs=0.01)
    assert set(rec["features"]) == set(names)


def test_parse_avi_response_missing_scorer():
    with pytest.raises(ValueError, match="lacks scorers"):
        scoring.parse_avi_response({"AVI_SCORE": fake_adata([1.0])})


def test_feature_columns_decomposition_and_loud_failure():
    names = all_feature_names()
    feats = {n: float(i) for i, n in enumerate(names)}
    am, ag = analysis.feature_columns(feats)
    assert am == feats[config.ALPHAMISSENSE_FEATURE]
    assert ag == sum(feats[n] for n in config.ALPHAGENOME_FEATURES)
    with pytest.raises(ValueError, match="not found"):
        analysis.feature_columns({"SOMETHING_ELSE": 1.0})
    # name formatting differences are tolerated
    pretty = {n.title().replace("_", " "): v for n, v in feats.items()}
    assert analysis.feature_columns(pretty)[0] == am


# ---- retry, budget, caching -------------------------------------------------- #
def test_retry_then_success_and_cache(tmp_path):
    cache = scoring.Cache(tmp_path)
    calls = {"n": 0}

    def flaky():
        calls["n"] += 1
        if calls["n"] < 3:
            raise TimeoutError("slow")
        return {"avi_raw": 1.0, "avi_quantile": 0.5, "avi_phred": 3.0, "features": {}}

    sleeps = []
    rec = scoring.score_variant_with_retry(cache, "chr1:1:A>G", flaky, sleep=sleeps.append)
    assert rec["avi_raw"] == 1.0 and calls["n"] == 3 and len(sleeps) == 2
    assert sleeps[1] > sleeps[0]  # exponential backoff
    assert cache.calls_made == 3
    # persisted: a new Cache object sees it and the counter
    again = scoring.Cache(tmp_path)
    assert "chr1:1:A>G" in again and again.calls_made == 3


def test_value_error_is_not_retried(tmp_path):
    cache = scoring.Cache(tmp_path)

    def bad():
        raise ValueError("invalid variant")

    assert scoring.score_variant_with_retry(cache, "k", bad, sleep=lambda s: None) is None
    assert cache.calls_made == 1


def test_permission_error_is_fatal(tmp_path):
    cache = scoring.Cache(tmp_path)

    def denied():
        raise PermissionError("secret-detail")

    with pytest.raises(scoring.FatalAPIError) as ei:
        scoring.score_variant_with_retry(cache, "k", denied, sleep=lambda s: None)
    assert "secret-detail" not in str(ei.value)


def test_budget_cap(tmp_path, monkeypatch):
    monkeypatch.setattr(config, "MAX_API_CALLS", 2)
    cache = scoring.Cache(tmp_path)
    cache.reserve_call()
    cache.reserve_call()
    with pytest.raises(scoring.BudgetExceeded):
        cache.reserve_call()


def test_api_key_never_in_repo_and_env_lookup(monkeypatch):
    monkeypatch.delenv("ALPHAGENOME_API_KEY", raising=False)
    monkeypatch.delenv("ALPHA_GENOME_API_KEY", raising=False)
    with pytest.raises(RuntimeError, match="No AlphaGenome API key"):
        scoring.get_api_key()
    monkeypatch.setenv("ALPHA_GENOME_API_KEY", "dummy-test-value")
    assert scoring.get_api_key() == "dummy-test-value"
