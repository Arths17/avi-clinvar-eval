#!/usr/bin/env bash
# Usage:  scripts/run_with_key.sh --dry-run          (20 variants)
#         scripts/run_with_key.sh --cohort all       (full run, after the dry run is checked)
# Asks for your AlphaGenome API key with hidden input, keeps it only in this process's
# environment, and runs `python -m avi_eval score` with the arguments you give.
set -euo pipefail
cd "$(dirname "$0")/.."
source .venv/bin/activate
read -r -s -p "Paste your NEW AlphaGenome API key (input is hidden), then press Enter: " ALPHAGENOME_API_KEY
echo
export ALPHAGENOME_API_KEY
python -m avi_eval score "$@"
