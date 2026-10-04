# AVI vs ClinVar: an independent evaluation

How well does the **AlphaGenome Variant Impact (AVI)** score from Google DeepMind's AlphaGenome Atlas separate
pathogenic from benign human variants in **ClinVar**, and how does that vary by consequence type (missense,
synonymous, splice, non-coding, ...)?

> **Research use only. No clinical or diagnostic claims. Nothing here is medical advice.**
> Independent researcher, no institutional affiliation. Not affiliated with or endorsed by Google or DeepMind.

The study plan was fixed *before* any score was retrieved: see [`design.md`](design.md).
All reported numbers come from code run on real data: aggregate results are in `results/` (`metrics.csv`, `summary.json`,
diagnostics, figures) and the paper in `paper/` is generated from those files by `paper/build_paper.py`
(`pip install -r requirements-paper.txt`, then `python paper/build_paper.py`). Per-variant AlphaGenome scores are not included.

## Repository layout

| Path | What |
|---|---|
| `design.md` | Pre-registered study design + dated deviations log |
| `src/avi_eval/` | `config.py` (all fixed parameters), `clinvar.py` (parse/filter), `sampling.py`, `scoring.py` (API + cache + retries), `metrics.py`, `analysis.py`, `plots.py`, `cli.py` |
| `tests/` | 55 unit tests (filters, labels, strata, sampling, metrics, bootstrap, retry/budget logic) on synthetic data |
| `notebooks/kaggle_avi_clinvar_eval.ipynb` | Kaggle-ready runner |
| `data/` | small derived files only (sampled ClinVar variant lists) |
| `results/`, `paper/` | outputs and write-up (filled in later stages) |
| `NOTICE_ALPHAGENOME_OUTPUT_TERMS.txt` | required AlphaGenome notice text |
| `requirements.txt` / `requirements-py311.txt` / `requirements-lock.txt` | pinned versions (Python 3.12 / 3.11 / full freeze) |

## Safety and terms (read this)

* **API key**: only from the environment variable `ALPHAGENOME_API_KEY` (or `ALPHA_GENOME_API_KEY`) or a Kaggle Secret.
  It is never printed, logged or written to disk. `.gitignore` excludes `.env`, `*.key`, `cache/` and `data/raw/`.
  Never paste your key into chat, a notebook cell, or a commit.
* **Per-variant AlphaGenome outputs** (`cache/`) stay local and out of git until the terms are re-checked (design §9).
  Keep any Kaggle notebook that has run the scoring **private**.
* **Non-commercial, no clinical use**, no training of variant-effect models on the outputs (the code only computes rank metrics).
* **Rate limits**: no numeric quota is published. The code uses 4 workers, exponential backoff, a disk cache, and a hard cap of
  6,000 API attempts (`config.MAX_API_CALLS`). The planned sample is at most ~4,860 variants.
* Notice for any published AlphaGenome Output: see `NOTICE_ALPHAGENOME_OUTPUT_TERMS.txt`.

## Run locally (Python 3.10+)

```bash
python3.12 -m venv .venv && source .venv/bin/activate
pip install -r requirements.txt && pip install -e . --no-deps
python -m pytest -q                        # 55 tests, no network needed

python -m avi_eval download                # shows file name and size; add --yes to download (~189 MB)
python -m avi_eval download --yes
python -m avi_eval prepare                 # parse + filter + sample; NO API calls
#   GATE: look at data/mc_term_counts.csv and results/attrition.csv first

export ALPHAGENOME_API_KEY="..."           # your personal key; do not share it
python -m avi_eval score --dry-run         # 20 variants end to end
python -m avi_eval score --cohort all      # full run (cached; safe to re-run)
python -m avi_eval analyze --n-boot 2000   # metrics, bootstrap CIs, figures
```

## Run on Kaggle: exactly what to click

1. Push this repository to your own GitHub (public is fine; **no keys, no `cache/`**).
2. On kaggle.com: **Create → New Notebook**. **File → Import Notebook** and upload `notebooks/kaggle_avi_clinvar_eval.ipynb`
   (or paste its cells).
3. Right sidebar → **Session options** → set **Internet** to **On** (Kaggle requires a phone-verified account for this).
4. Menu **Add-ons → Secrets → Add Secret**. Label: `ALPHAGENOME_API_KEY`. Value: your key. Then make sure the secret's
   toggle/checkbox is **on for this notebook** ("Attached").
5. Make sure the notebook is **Private** (Share button → Private).
6. In the second cell set `REPO_URL` to your GitHub repo URL. Leave `DRY_RUN = True`.
7. **Run → Run All.** Read the printed attrition table, the MC-term table and the 20-variant dry-run table.
   Send those to Claude. Do **not** continue until the gate in `design.md` is checked.
8. Change `DRY_RUN = False` and **Run All** again (the full run reuses the disk cache and takes on the order of tens of minutes).
9. Download `results_bundle.zip` from the notebook's **Output** tab and share `results/metrics.csv` with Claude.

## Data sources and citations

* ClinVar (NCBI): Landrum et al., *Nucleic Acids Res* (see NCBI's ClinVar citation page); file used: `clinvar_20260928.vcf.gz`.
* AlphaGenome: Avsec et al., "Advancing regulatory variant effect prediction with AlphaGenome", *Nature* (2026), doi:10.1038/s41586-025-10014-0.
* AlphaGenome Atlas / AVI: Cheng et al., "AlphaGenome Atlas: in silico mutagenesis of the entire human genome improves prioritization and
  interpretation of non-coding variants", medRxiv (preprint), doi:10.64898/2026.09.16.26363192.

## AI assistance disclosure

Code and drafting were assisted by Claude (Anthropic). The author reviews and takes responsibility for the content.

## License

Code: MIT (see `LICENSE`). AlphaGenome outputs and ClinVar data are not covered by this license.
