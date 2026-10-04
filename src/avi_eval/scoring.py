"""Retrieve AVI scores from the AlphaGenome Atlas API, safely.

Safety properties (design.md section 8):
  * The API key is read from an environment variable or Kaggle Secrets. It is never
    printed, logged, written to disk or put in an exception message.
  * Every successful response is appended to a JSONL cache on disk; reruns make no
    new calls for cached variants.
  * A persistent call counter enforces a hard cap (config.MAX_API_CALLS).
  * Transient errors are retried with exponential backoff and jitter.
"""
from __future__ import annotations

import json
import logging
import math
import os
import random
import threading
import time
from concurrent.futures import ThreadPoolExecutor, as_completed
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Callable, Mapping

import numpy as np
import pandas as pd

from . import config

log = logging.getLogger(__name__)


class BudgetExceeded(RuntimeError):
    pass


class FatalAPIError(RuntimeError):
    """Authentication / permission / not-implemented errors: retrying cannot help."""


# --------------------------------------------------------------------------- #
# Keys and cache keys
# --------------------------------------------------------------------------- #
def variant_key(chrom: str, pos: int, ref: str, alt: str) -> str:
    return f"chr{chrom}:{int(pos)}:{ref}>{alt}"


def get_api_key() -> str:
    """Return the key from env vars, falling back to Kaggle Secrets. Never prints it."""
    for name in config.API_KEY_ENV_NAMES:
        val = os.environ.get(name)
        if val:
            return val
    try:  # Kaggle notebooks
        from kaggle_secrets import UserSecretsClient  # type: ignore

        client = UserSecretsClient()
        for name in config.API_KEY_ENV_NAMES:
            try:
                val = client.get_secret(name)
            except Exception:  # secret not defined under this name
                continue
            if val:
                return val
    except ImportError:
        pass
    raise RuntimeError(
        "No AlphaGenome API key found. Set the environment variable "
        f"{config.API_KEY_ENV_NAMES[0]} (or add a Kaggle Secret with that name). "
        "Never paste the key into code or chat."
    )


# --------------------------------------------------------------------------- #
# Disk cache + call budget
# --------------------------------------------------------------------------- #
class Cache:
    def __init__(self, cache_dir: Path = config.CACHE_DIR):
        self.dir = Path(cache_dir)
        self.dir.mkdir(parents=True, exist_ok=True)
        self.path = self.dir / "avi_cache.jsonl"
        self.error_path = self.dir / "avi_errors.jsonl"
        self.counter_path = self.dir / "api_calls.json"
        self._lock = threading.Lock()
        self._records: dict[str, dict] = {}
        if self.path.exists():
            with open(self.path, encoding="utf-8") as fh:
                for line in fh:
                    line = line.strip()
                    if line:
                        rec = json.loads(line)
                        self._records[rec["key"]] = rec
        self._calls = 0
        if self.counter_path.exists():
            self._calls = int(json.loads(self.counter_path.read_text())["calls"])

    def __contains__(self, key: str) -> bool:
        return key in self._records

    def get(self, key: str) -> dict | None:
        return self._records.get(key)

    def records(self) -> dict[str, dict]:
        return dict(self._records)

    @property
    def calls_made(self) -> int:
        return self._calls

    def reserve_call(self) -> None:
        """Count one API attempt before making it; raise if the hard cap would be passed."""
        with self._lock:
            if self._calls >= config.MAX_API_CALLS:
                raise BudgetExceeded(
                    f"Hard API-call cap of {config.MAX_API_CALLS} reached (design.md section 4). Stopping."
                )
            self._calls += 1
            self.counter_path.write_text(json.dumps({"calls": self._calls}))

    def put(self, rec: dict) -> None:
        with self._lock:
            self._records[rec["key"]] = rec
            with open(self.path, "a", encoding="utf-8") as fh:
                fh.write(json.dumps(rec) + "\n")

    def put_error(self, key: str, message: str) -> None:
        with self._lock:
            with open(self.error_path, "a", encoding="utf-8") as fh:
                fh.write(json.dumps({"key": key, "error": message[:300], "time": _now()}) + "\n")


def _now() -> str:
    return datetime.now(timezone.utc).isoformat(timespec="seconds")


# --------------------------------------------------------------------------- #
# Parsing the API response
# --------------------------------------------------------------------------- #
def phred_from_cdf_quantile(q: float) -> float:
    """PHRED = -10*log10(1 - cdf_quantile), the convention in DeepMind's science-skills AVI guide.

    VERIFY on the dry run against the Atlas website for 1-2 variants; the raw score and the
    quantile are stored so PHRED can be recomputed if this convention turns out to differ.
    """
    tail = max(1.0 - float(q), 1e-12)
    return float(-10.0 * math.log10(tail))


def _feature_names(adata) -> list[str]:
    var = adata.var
    if "name" in var.columns:
        return [str(x) for x in var["name"]]
    return [str(x) for x in adata.var_names]


def parse_avi_response(scores: Mapping[str, Any]) -> dict:
    """Turn the {scorer_name: AnnData} mapping into a flat, JSON-serialisable record."""
    missing = [s for s in config.REQUESTED_SCORERS if s not in scores]
    if missing:
        raise ValueError(f"Response lacks scorers {missing}; got {sorted(scores)}")
    avi = scores["AVI_SCORE"]
    x = np.ravel(np.asarray(avi.X, dtype=float))
    if x.size != 1:
        raise ValueError(f"Expected exactly one AVI_SCORE value, got shape {avi.X.shape}")
    quantile = None
    if "quantiles" in avi.layers:
        quantile = float(np.ravel(np.asarray(avi.layers["quantiles"], dtype=float))[0])
    fi = scores["AVI_SCORE_FEATURE_IMPORTANCE"]
    values = np.ravel(np.asarray(fi.X, dtype=float))
    names = _feature_names(fi)
    if len(names) != len(values):
        raise ValueError(f"Feature importance has {len(values)} values but {len(names)} names")
    return {
        "avi_raw": float(x[0]),
        "avi_quantile": quantile,
        "avi_phred": None if quantile is None else phred_from_cdf_quantile(quantile),
        "features": {n: float(v) for n, v in zip(names, values)},
    }


# --------------------------------------------------------------------------- #
# Calling the API
# --------------------------------------------------------------------------- #
def make_client():  # pragma: no cover - needs network and a key
    from alphagenome.atlas import atlas

    return atlas.create(get_api_key())


def query_one(client, chrom: str, pos: int, ref: str, alt: str) -> dict:  # pragma: no cover - network
    from alphagenome.data import genome

    variant = genome.Variant(chromosome=f"chr{chrom}", position=int(pos), reference_bases=ref, alternate_bases=alt)
    scores = client.query_variant(variant, requested_scorers=list(config.REQUESTED_SCORERS))
    return parse_avi_response(scores)


def _classify_error(exc: BaseException) -> str:
    """'fatal' | 'retry' | 'skip'."""
    if isinstance(exc, (PermissionError, NotImplementedError)):
        return "fatal"
    if isinstance(exc, (ValueError, IndexError)):
        return "skip"  # bad request for this variant; retrying will not help
    return "retry"  # timeouts, gRPC unavailable / resource exhausted, etc.


def score_variant_with_retry(
    cache: Cache,
    key: str,
    call: Callable[[], dict],
    *,
    sleep: Callable[[float], None] = time.sleep,
    max_attempts: int = config.MAX_ATTEMPTS,
) -> dict | None:
    """Call the API for one variant with budget accounting and backoff. Returns the record or None."""
    for attempt in range(1, max_attempts + 1):
        cache.reserve_call()
        try:
            result = call()
        except BudgetExceeded:
            raise
        except Exception as exc:  # noqa: BLE001 - classified below
            kind = _classify_error(exc)
            cache.put_error(key, f"{type(exc).__name__}: {exc}")
            if kind == "fatal":
                raise FatalAPIError(f"{type(exc).__name__} for {key} - check key/permissions") from None
            if kind == "skip" or attempt == max_attempts:
                log.warning("giving up on %s after %d attempt(s): %s", key, attempt, type(exc).__name__)
                return None
            delay = config.BACKOFF_BASE_SECONDS * (2 ** (attempt - 1)) + random.uniform(0, 1)
            log.info("retrying %s in %.1fs (%s)", key, delay, type(exc).__name__)
            sleep(delay)
            continue
        rec = {"key": key, "time": _now(), "alphagenome_pkg": _pkg_version(), **result}
        cache.put(rec)
        return rec
    return None


def _pkg_version() -> str:
    try:
        from importlib.metadata import version

        return version("alphagenome")
    except Exception:  # noqa: BLE001
        return "unknown"


def score_dataframe(
    sample: pd.DataFrame,
    cache: Cache,
    *,
    client=None,
    workers: int = config.WORKERS,
) -> dict:  # pragma: no cover - network orchestration
    """Score every variant in `sample` that is not already cached. Returns a small summary."""
    keys = [variant_key(r.chrom, r.pos, r.ref, r.alt) for r in sample.itertuples()]
    todo = [(k, r) for k, r in zip(keys, sample.itertuples()) if k not in cache]
    log.info("%d variants in sample, %d already cached, %d to query", len(keys), len(keys) - len(todo), len(todo))
    if not todo:
        return {"queried": 0, "ok": 0, "failed": 0}
    client = client or make_client()
    ok = failed = 0
    with ThreadPoolExecutor(max_workers=workers) as pool:
        futs = {
            pool.submit(
                score_variant_with_retry, cache, k,
                lambda r=r: query_one(client, r.chrom, r.pos, r.ref, r.alt),
            ): k
            for k, r in todo
        }
        try:
            from tqdm import tqdm

            it = tqdm(as_completed(futs), total=len(futs), desc="AVI queries")
        except ImportError:
            it = as_completed(futs)
        for fut in it:
            try:
                rec = fut.result()
            except (BudgetExceeded, FatalAPIError):
                for f in futs:
                    f.cancel()
                raise
            if rec is None:
                failed += 1
            else:
                ok += 1
    return {"queried": len(todo), "ok": ok, "failed": failed}


def attach_scores(sample: pd.DataFrame, cache: Cache) -> pd.DataFrame:
    """Join cached AVI results onto the sample (rows without a cached score are kept with NaN)."""
    recs = cache.records()
    rows = []
    for r in sample.itertuples():
        rec = recs.get(variant_key(r.chrom, r.pos, r.ref, r.alt))
        rows.append(
            {
                "avi_raw": np.nan if rec is None else rec["avi_raw"],
                "avi_quantile": np.nan if rec is None or rec["avi_quantile"] is None else rec["avi_quantile"],
                "avi_phred": np.nan if rec is None or rec["avi_phred"] is None else rec["avi_phred"],
                "features": None if rec is None else rec["features"],
            }
        )
    return pd.concat([sample.reset_index(drop=True), pd.DataFrame(rows)], axis=1)
