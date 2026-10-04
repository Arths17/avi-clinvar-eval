"""Command line: download -> prepare -> score (dry run first!) -> analyze."""
from __future__ import annotations

import argparse
import hashlib
import json
import logging
import platform
import sys
import urllib.request
from datetime import datetime, timezone
from importlib.metadata import PackageNotFoundError, version
from pathlib import Path

import pandas as pd

from . import analysis, clinvar, config, plots, sampling, scoring

log = logging.getLogger("avi_eval")


def _setup_logging() -> None:
    config.RESULTS_DIR.mkdir(parents=True, exist_ok=True)
    handlers = [logging.StreamHandler(sys.stdout), logging.FileHandler(config.RESULTS_DIR / "run.log")]
    logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s", handlers=handlers)


def _pkg(name: str) -> str:
    try:
        return version(name)
    except PackageNotFoundError:
        return "not installed"


def _hash_file(path: Path) -> tuple[str, str]:
    md5, sha = hashlib.md5(), hashlib.sha256()  # noqa: S324 - md5 only to match NCBI's checksum file
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            md5.update(chunk)
            sha.update(chunk)
    return md5.hexdigest(), sha.hexdigest()


def _vcf_path(release: str) -> Path:
    return config.RAW_DIR / f"clinvar_{release}.vcf.gz"


def _update_manifest(update: dict) -> None:
    path = config.RESULTS_DIR / "run_manifest.json"
    cur = json.loads(path.read_text()) if path.exists() else {}
    cur.update(update)
    path.write_text(json.dumps(cur, indent=2, default=str))


# --------------------------------------------------------------------------- #
def cmd_download(args) -> None:
    rel = args.release
    url = f"{config.CLINVAR_BASE_URL}/clinvar_{rel}.vcf.gz"
    dest = _vcf_path(rel)
    config.RAW_DIR.mkdir(parents=True, exist_ok=True)
    with urllib.request.urlopen(url) as r:  # noqa: S310 - fixed NCBI https URL
        size_mb = int(r.headers.get("Content-Length", 0)) / 1e6
    print(f"Will download {url}\n  to {dest}\n  size ~{size_mb:.0f} MB (ClinVar, NCBI; public data).")
    if not args.yes:
        print("Re-run with --yes to confirm the download.")
        return
    with urllib.request.urlopen(url) as r, open(dest, "wb") as out:  # noqa: S310
        for chunk in clinvar.iter_download_chunks(r):
            out.write(chunk)
    with urllib.request.urlopen(url + ".md5") as r:  # noqa: S310
        expected = r.read().decode().split()[0].strip().lower()
    md5, sha = _hash_file(dest)
    if md5 != expected:
        dest.unlink()
        raise SystemExit(f"MD5 mismatch (got {md5}, NCBI says {expected}); file deleted.")
    print(f"MD5 verified: {md5}")
    _update_manifest({"clinvar_file": dest.name, "clinvar_md5": md5, "clinvar_sha256": sha, "clinvar_url": url})


def cmd_prepare(args) -> None:
    vcf = Path(args.vcf) if args.vcf else _vcf_path(config.CLINVAR_RELEASE)
    att = clinvar.Attrition()
    df, unknown = clinvar.load_clinvar(vcf, att)
    if unknown:
        log.warning("Unrecognised review-status strings (excluded): %s", dict(unknown.most_common(10)))
    df = clinvar.add_strata(df)
    mc = clinvar.mc_term_counts(df)
    config.DATA_DIR.mkdir(parents=True, exist_ok=True)
    mc.to_csv(config.DATA_DIR / "mc_term_counts.csv", index=False)
    unmapped = mc[mc["mapped_to"] == "UNMAPPED"]
    if len(unmapped):
        log.warning("MC terms not mapped to any stratum (see data/mc_term_counts.csv):\n%s", unmapped.head(15).to_string(index=False))
    primary, onestar = sampling.build_primary_and_onestar(df, att)
    primary.to_csv(config.DATA_DIR / "sample_primary.csv", index=False)
    onestar.to_csv(config.DATA_DIR / "sample_onestar.csv", index=False)
    sampling.pick_dry_run(primary).to_csv(config.DATA_DIR / "sample_dry_run.csv", index=False)
    att.save(config.RESULTS_DIR / "attrition.csv")
    total = len(primary) + len(onestar)
    print(f"\nPrimary sample: {len(primary)} variants; 1-star supplement: {len(onestar)}; total {total}.")
    print("GATE: review data/mc_term_counts.csv and results/attrition.csv BEFORE running `score`.")
    _update_manifest(
        {
            "prepared_at": datetime.now(timezone.utc).isoformat(timespec="seconds"),
            "seed": config.SEED, "python": platform.python_version(),
            "pandas": _pkg("pandas"), "numpy": _pkg("numpy"), "scikit-learn": _pkg("scikit-learn"),
            "alphagenome": _pkg("alphagenome"), "n_primary": len(primary), "n_onestar": len(onestar),
        }
    )


def cmd_score(args) -> None:
    cache = scoring.Cache(config.CACHE_DIR)
    if args.dry_run:
        sample = pd.read_csv(config.DATA_DIR / "sample_dry_run.csv")
    else:
        files = {"primary": ["sample_primary.csv"], "onestar": ["sample_onestar.csv"],
                 "all": ["sample_primary.csv", "sample_onestar.csv"]}[args.cohort]
        sample = pd.concat([pd.read_csv(config.DATA_DIR / f) for f in files], ignore_index=True)
        sample = sample.drop_duplicates("variation_id")
    pending = sum(scoring.variant_key(r.chrom, r.pos, r.ref, r.alt) not in cache for r in sample.itertuples())
    print(f"{len(sample)} variants selected; {pending} not yet cached; API calls used so far: {cache.calls_made}/{config.MAX_API_CALLS}.")
    if pending + cache.calls_made > config.MAX_API_CALLS:
        raise SystemExit("This run could exceed the hard API-call cap; refusing to start.")
    summary = scoring.score_dataframe(sample, cache)
    print("Scoring summary:", summary)
    scored = scoring.attach_scores(sample, cache)
    if args.dry_run:
        out = config.RESULTS_DIR / "dry_run"
        out.mkdir(parents=True, exist_ok=True)
        show = scored.drop(columns=["features"])
        show.to_csv(out / "dry_run_scores.csv", index=False)
        print(show[["variation_id", "chrom", "pos", "ref", "alt", "label", "stratum", "avi_raw", "avi_quantile", "avi_phred"]].to_string(index=False))
        feats = next((f for f in scored["features"] if f), None)
        if feats:
            print("\nFeature names returned by the API (check against config.ALPHAGENOME_FEATURES):")
            print(sorted(feats))
        print("\nNext: compare AVI/PHRED for 1-2 variants on the Atlas website (single-variant lookup),"
              " and check that no variant failed.")
    _update_manifest({"api_calls_used": cache.calls_made})


def cmd_analyze(args) -> None:
    cache = scoring.Cache(config.CACHE_DIR)
    onestar = config.DATA_DIR / "sample_onestar.csv"
    res = analysis.run_all(config.DATA_DIR / "sample_primary.csv", onestar, cache, config.RESULTS_DIR, args.n_boot)
    plots.make_all(res, config.RESULTS_DIR / "figures")
    s = res["summary"]
    print("\nSanity gate:", "INVESTIGATE BEFORE INTERPRETING" if s["sanity"]["needs_investigation"] else "no flags")
    for f in s["sanity"]["flags"]:
        print("  -", f)
    print("Outputs written to results/ (metrics.csv, summary.json, figures/). Paste results/metrics.csv back to Claude.")


def main(argv: list[str] | None = None) -> None:
    p = argparse.ArgumentParser(prog="avi_eval", description=__doc__)
    sub = p.add_subparsers(dest="cmd", required=True)

    d = sub.add_parser("download", help="download the ClinVar VCF (asks for --yes)")
    d.add_argument("--release", default=config.CLINVAR_RELEASE)
    d.add_argument("--yes", action="store_true")
    d.set_defaults(fn=cmd_download)

    pr = sub.add_parser("prepare", help="parse ClinVar, filter, sample (no API calls)")
    pr.add_argument("--vcf", default=None)
    pr.set_defaults(fn=cmd_prepare)

    sc = sub.add_parser("score", help="query AVI from the Atlas API (cached)")
    sc.add_argument("--dry-run", action="store_true", help=f"only {config.DRY_RUN_N} variants")
    sc.add_argument("--cohort", choices=["primary", "onestar", "all"], default="primary")
    sc.set_defaults(fn=cmd_score)

    an = sub.add_parser("analyze", help="compute metrics, bootstrap CIs and figures from cached scores")
    an.add_argument("--n-boot", type=int, default=config.N_BOOT)
    an.set_defaults(fn=cmd_analyze)

    args = p.parse_args(argv)
    _setup_logging()
    args.fn(args)


if __name__ == "__main__":
    main()
