#!/usr/bin/env bash
set -euo pipefail

bash scripts/install_tape_cpu.sh
bash scripts/install_blade_py311.sh
bash scripts/install_feature_selector_py311.sh
bash scripts/install_blue_cpu.sh
