#!/usr/bin/env bash
set -euo pipefail

CHANNEL_ARGS=(
  -c conda-forge
  -c bioconda
)

DRY_RUN_ENV_NAME=foli-decon-expansion-dry

PACKAGE_CHECKS=(
  "TAPE_R_ALIAS:r-tape"
  "BLADE_R_ALIAS:r-blade"
  "MEAD:r-mead"
  "ABIS:r-abis"
  "ESTIMATE:r-estimate"
  "TIMER:r-timer"
  "TIMER2:r-timer2"
  "TIMER3:r-timer3"
  "CELLCODE:r-cellcode"
  "DSA:r-dsa"
  "BLUE:r-blue"
  "PSEA:bioconductor-psea=1.32.0"
  "AUTOGENES:autogenes=1.0.4"
  "SCGENEFIT:scgenefit"
  "MARKERMAP:markermap"
)

PIP_CHECKS=(
  "TAPE:scTAPE:scTAPE==1.1.2"
  "BLADE:BLADE-Deconvolution:BLADE-Deconvolution==0.0.7"
  "SCGENEFIT:scGeneFit:scGeneFit==1.0.2"
  "MARKERMAP:markermap:markermap==1.0.3"
)

PIP_FALSE_POSITIVE_CHECKS=(
  "MEAD:mead"
  "BLUE:blue"
  "DSA:dsa"
  "TIMER:timer"
)

SIDECAR_ENV_CHECKS=(
  "R44_MERGED:environment-r44.yml"
  "PY311_CPU_MERGED:environment-py311-cpu.yml"
  "SCDC_R40:environment-scdc-r40.yml"
  "DEEP_CPU:environment-deep-cpu.yml"
)

echo "=== Requested method package search (conda-forge + bioconda) ==="
for entry in "${PACKAGE_CHECKS[@]}"; do
  name="${entry%%:*}"
  spec="${entry#*:}"
  echo "${name} :: mamba search -c conda-forge -c bioconda ${spec}"
  if mamba search "${CHANNEL_ARGS[@]}" "${spec}" | grep -qi "No entries matching\|No match"; then
    echo "  result: not found in channel index"
  else
    echo "  result: present"
  fi
  echo
done

echo "=== Requested method source/package probes outside conda ==="
echo "ESTIMATE :: R-Forge estimate package"
Rscript -e 'rforge <- "http://r-forge.r-project.org"; ap <- available.packages(repos=rforge); hit <- intersect("estimate", rownames(ap)); if (length(hit) > 0) { print(ap[hit, intersect(c("Package", "Version", "Depends", "Imports", "NeedsCompilation"), colnames(ap)), drop=FALSE], quote=FALSE) } else { cat("  result: not found on R-Forge\n") }'
echo
echo "ABIS :: GitHub Shiny source resources"
python - <<'PY'
from urllib.request import urlopen

for path in ("server.R", "data/sigmatrixRNAseq.txt", "data/sigmatrixMicro.txt", "data/target.txt"):
    url = f"https://raw.githubusercontent.com/giannimonaco/ABIS/master/{path}"
    with urlopen(url, timeout=30) as response:
        print(f"  {path}: {response.headers.get('Content-Length', 'unknown')} bytes")
PY
echo
echo "MEAD :: GitHub R package source route"
python - <<'PY'
from urllib.request import urlopen

for url in (
    "https://raw.githubusercontent.com/DongyueXie/MEAD/master/DESCRIPTION",
    "https://raw.githubusercontent.com/mengyin/vashr/master/DESCRIPTION",
):
    with urlopen(url, timeout=30) as response:
        first_line = response.readline().decode("utf-8").strip()
        print(f"  {url}: {first_line}")
PY
echo
echo "CellCODE :: GitHub R package source route"
python - <<'PY'
from urllib.request import urlopen

url = "https://raw.githubusercontent.com/mchikina/CellCODE/master/DESCRIPTION"
with urlopen(url, timeout=30) as response:
    first_line = response.readline().decode("utf-8").strip()
    print(f"  {url}: {first_line}")
PY
echo
echo "BLUE :: GitHub Python pipeline source route"
python - <<'PY'
from urllib.request import urlopen

for path in ("pyproject.toml", "README.md", "CustomDS_README.md"):
    url = f"https://raw.githubusercontent.com/SichenZhu/BLUE/main/{path}"
    with urlopen(url, timeout=30) as response:
        first_line = response.readline().decode("utf-8").strip()
        print(f"  {path}: {first_line}")
PY
echo

echo "=== Requested method package search (PyPI) ==="
for entry in "${PIP_CHECKS[@]}"; do
  name="${entry%%:*}"
  rest="${entry#*:}"
  package="${rest%%:*}"
  echo "${name} :: python -m pip index versions ${package}"
  if python -m pip index versions "$package" >/tmp/requested-method-pip-index.log 2>&1; then
    sed -n '1,3p' /tmp/requested-method-pip-index.log
  else
    echo "  result: not found on PyPI"
  fi
  echo
done

echo "=== Requested unresolved-name PyPI false-positive checks ==="
for entry in "${PIP_FALSE_POSITIVE_CHECKS[@]}"; do
  name="${entry%%:*}"
  package="${entry#*:}"
  echo "${name} :: python -m pip install --dry-run ${package}"
  if python -m pip install --dry-run "$package" >/tmp/requested-method-pip-false-positive.log 2>&1; then
    grep -E "Would install|Collecting" /tmp/requested-method-pip-false-positive.log | sed -n '1,20p'
  else
    echo "  result: not found on PyPI"
  fi
  echo
done

echo "=== Requested method dry-run install probes (R 4.2 where possible) ==="
DRY_RUN_LIST=(
  "r-base=4.2 r-tape"
  "r-base=4.2 r-blade"
  "r-base=4.2 r-mead"
  "r-base=4.2 r-abis"
  "r-base=4.2 r-estimate"
  "r-base=4.2 r-mass bioconductor-preprocesscore"
  "r-base=4.2 r-cellcode"
  "r-base=4.2 r-dsa"
  "r-base=4.2 r-blue"
  "r-base=4.2 bioconductor-psea=1.32.0"
  "python=3.11 autogenes=1.0.4"
  "python=3.11 scgenefit"
)
for spec in "${DRY_RUN_LIST[@]}"; do
  echo "  mamba create --dry-run -n ${DRY_RUN_ENV_NAME} ${CHANNEL_ARGS[*]} ${spec}"
  if mamba create --dry-run -n "$DRY_RUN_ENV_NAME" "${CHANNEL_ARGS[@]}" $spec >/tmp/requested-method-dryrun.log 2>&1; then
    echo "  result: solved"
  else
    echo "  result: failed to solve"
    grep -E "package.*could not be installed|Could not set|Could not solve|No entries" /tmp/requested-method-dryrun.log || true
  fi
  echo
 done

echo "=== Requested method dry-run install probes (PyPI only; no install) ==="
for entry in "${PIP_CHECKS[@]}"; do
  name="${entry%%:*}"
  rest="${entry#*:}"
  spec="${rest#*:}"
  echo "  python -m pip install --dry-run ${spec}"
  if python -m pip install --dry-run "$spec" >/tmp/requested-method-pip-dryrun.log 2>&1; then
    grep -E "Would install|Requirement already satisfied|Successfully downloaded|Using cached|Collecting" /tmp/requested-method-pip-dryrun.log | sed -n '1,30p'
  else
    echo "  result: failed"
    sed -n '1,80p' /tmp/requested-method-pip-dryrun.log
  fi
  echo
done

echo "=== Requested method sidecar env dry-runs ==="
for entry in "${SIDECAR_ENV_CHECKS[@]}"; do
  name="${entry%%:*}"
  file="${entry#*:}"
  echo "  ${name} :: mamba env create --dry-run -f ${file}"
  if mamba env create --dry-run -f "$file" >/tmp/requested-method-sidecar-dryrun.log 2>&1; then
    echo "  result: solved"
  else
    echo "  result: failed"
    grep -E "package.*could not be installed|Could not set|Could not solve|No entries" /tmp/requested-method-sidecar-dryrun.log || true
  fi
  echo
done

echo "=== Sweep complete ==="
