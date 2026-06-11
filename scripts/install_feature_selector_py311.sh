#!/usr/bin/env bash
set -euo pipefail

python -m pip install --no-deps scGeneFit==1.0.2
python -m pip install --no-deps lassonet==0.0.20
python -m pip install --no-deps git+https://github.com/iancovert/persist
python -m pip install --no-deps markermap==1.0.3
