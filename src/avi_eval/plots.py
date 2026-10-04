"""Static figures for the paper. Okabe-Ito colour-blind-safe colours, one axis per chart, thin marks."""
from __future__ import annotations

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from sklearn.metrics import roc_curve  # noqa: E402

from . import config  # noqa: E402

BLUE, ORANGE, GREY = "#0072B2", "#E69F00", "#6b6b6b"
INK = "#1a1a1a"

plt.rcParams.update(
    {
        "font.size": 10, "axes.edgecolor": GREY, "axes.labelcolor": INK, "text.color": INK,
        "xtick.color": INK, "ytick.color": INK, "axes.spines.top": False, "axes.spines.right": False,
        "axes.grid": True, "grid.color": "#e3e3e3", "grid.linewidth": 0.6, "figure.dpi": 150,
        "savefig.bbox": "tight",
    }
)


def roc_pooled(df: pd.DataFrame, path: Path) -> None:
    fig, ax = plt.subplots(figsize=(4.6, 4.4))
    for col, name, color in (("avi_raw", "AVI", BLUE), ("severity_b1", "Consequence baseline (B1)", ORANGE)):
        fpr, tpr, _ = roc_curve(df["label"], df[col])
        ax.plot(fpr, tpr, color=color, lw=2, label=name)
    ax.plot([0, 1], [0, 1], color=GREY, lw=1, ls="--", label="Chance")
    ax.set_xlabel("False-positive rate (benign called pathogenic)")
    ax.set_ylabel("True-positive rate (pathogenic found)")
    ax.set_title("Pooled ROC curve, primary sample", loc="left", fontsize=10)
    ax.legend(frameon=False, loc="lower right")
    ax.set_aspect("equal")
    fig.savefig(path)
    plt.close(fig)


def auc_by_stratum(table: pd.DataFrame, path: Path, analysis_set: str = "primary") -> None:
    t = table[(table["analysis_set"] == analysis_set) & (table["score"] == "avi_raw")
              & table["group"].str.startswith("stratum:") & table["eligible_for_auc"]].copy()
    if t.empty:
        return
    t["name"] = t["group"].str.replace("stratum:", "", regex=False)
    t = t.sort_values("roc_auc")
    fig, ax = plt.subplots(figsize=(6.0, 0.55 * len(t) + 1.4))
    y = np.arange(len(t))
    ax.hlines(y, t["roc_lo"], t["roc_hi"], color=BLUE, lw=1.5)
    ax.plot(t["roc_auc"], y, "o", color=BLUE, ms=7, mec="white", mew=1.5)
    ax.axvline(0.5, color=GREY, lw=1, ls="--")
    ax.set_yticks(y)
    ax.set_yticklabels([f"{n}  (n={int(p)}+/{int(q)}-)" for n, p, q in zip(t["name"], t["n_pos"], t["n_neg"])])
    ax.set_xlim(0.4, 1.0)
    ax.set_xlabel("ROC-AUC (95% gene-clustered bootstrap CI); dashed = chance")
    ax.set_title("AVI discrimination by consequence category", loc="left", fontsize=10)
    ax.grid(axis="y", visible=False)
    fig.savefig(path)
    plt.close(fig)


def class_counts(counts: pd.DataFrame, path: Path) -> None:
    fig, ax = plt.subplots(figsize=(6.0, 3.6))
    x = np.arange(len(counts))
    w = 0.38
    ax.bar(x - w / 2, counts["n_pathogenic"], w - 0.04, color=ORANGE, label="Pathogenic / likely pathogenic")
    ax.bar(x + w / 2, counts["n_benign"], w - 0.04, color=BLUE, label="Benign / likely benign")
    ax.set_xticks(x)
    ax.set_xticklabels(counts.index, rotation=35, ha="right")
    ax.set_ylabel("Variants in sample")
    ax.set_title("Class balance by consequence category (primary sample)", loc="left", fontsize=10)
    ax.legend(frameon=False)
    ax.grid(axis="x", visible=False)
    fig.savefig(path)
    plt.close(fig)


def score_distribution(df: pd.DataFrame, path: Path) -> None:
    fig, ax = plt.subplots(figsize=(5.4, 3.6))
    bins = np.linspace(df["avi_phred"].min(), df["avi_phred"].max(), 30)
    ax.hist(df.loc[df["label"] == 0, "avi_phred"], bins=bins, alpha=0.6, color=BLUE, label="Benign / likely benign")
    ax.hist(df.loc[df["label"] == 1, "avi_phred"], bins=bins, alpha=0.6, color=ORANGE, label="Pathogenic / likely pathogenic")
    for t in config.PHRED_THRESHOLDS:
        ax.axvline(t, color=GREY, lw=0.8, ls=":")
    ax.set_xlabel("AVI PHRED (dotted lines: 10, 20, 30)")
    ax.set_ylabel("Variants")
    ax.legend(frameon=False)
    ax.set_title("AVI PHRED by class, primary sample", loc="left", fontsize=10)
    fig.savefig(path)
    plt.close(fig)


def make_all(result: dict, out_dir: Path) -> None:
    out_dir.mkdir(parents=True, exist_ok=True)
    roc_pooled(result["primary"], out_dir / "roc_pooled.png")
    auc_by_stratum(result["table"], out_dir / "auc_by_stratum.png")
    class_counts(result["counts"], out_dir / "class_counts.png")
    score_distribution(result["primary"], out_dir / "avi_phred_distribution.png")
