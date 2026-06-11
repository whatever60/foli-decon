#!/usr/bin/env bash
set -euo pipefail

BLUE_REPO="${FOLI_DECON_BLUE_REPO:-resources/blue/BLUE}"

if [ ! -d "${BLUE_REPO}/.git" ]; then
  mkdir -p "$(dirname "${BLUE_REPO}")"
  git clone --depth 1 https://github.com/SichenZhu/BLUE.git "${BLUE_REPO}"
fi

python - <<'PY'
import anndata
import pandas
import scanpy
import torch
import yaml

print("BLUE CPU sidecar imports ok")
print(f"torch={torch.__version__}")
PY
