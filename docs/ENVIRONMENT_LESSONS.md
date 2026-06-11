# Environment Creation Lessons

This note records the dependency lessons from building the unified `foli-decon`
environment. The short version is that the environment is viable as one conda
env, but only if the heavy R stack is installed from conda first and source
installs are kept very narrow.

## Current Working Strategy

Use this order:

1. Create `foli-decon` from `environment.yml`.
2. Let conda/mamba install all heavy R packages and compiled dependencies.
3. Run `scripts/install_r_packages.R` only for the remaining R packages that are
   not available in the conda solve.

The main env and every sidecar YAML were dry-run solved again on 2026-06-11:

- `environment.yml`
- `environment-r44.yml`
- `environment-scdc-r40.yml`
- `environment-deep-cpu.yml`
- `environment-py311-cpu.yml`

The compressed sidecar layout is:

- `foli-decon-r44`: MuSiC 1.0.0/MuSiC2, MEAD, CellCODE, InstaPrism, CDSeq
- `foli-decon-scdc-r40`: SCDC only
- `foli-decon-deep-cpu`: Scaden and DISSECT
- `foli-decon-py311-cpu`: AutoGeneS, BLADE, BLUE, TAPE, MarkerMap, scGeneFit

Older per-tool sidecar YAMLs are retained as troubleshooting/probe artifacts,
but they are not the preferred installation route.

The local probe logs live under `cache/mamba_probe/` and are intentionally
ignored by git because they are machine-specific solver artifacts.

The post-install script intentionally checks that the conda-managed R packages
already exist before doing anything else. That guard prevents accidentally
triggering a large recursive CRAN/Bioconductor source build.

## Main Hurdles

### BayesPrism controls the R version

`r-bayesprism` is the package that most strongly constrains the unified env.
The working env uses `r-base=4.2` because the conda `r-bayesprism` build is tied
to the R 4.2 stack. Trying to push the environment to a newer R should be treated
as a full environment re-solve, not a small upgrade.

Current working snapshot:

- R: `4.2.3`
- conda package: `r-bayesprism`
- R package loaded as: `BayesPrism 2.0`

Practical lesson: prioritize BayesPrism compatibility over "latest R" if we want
all methods in one environment.

### MuSiC exact replication needs an R 4.4 sidecar

The NEC MuSiC notebook was originally run in the separate `folitools` env, not
in the unified `foli-decon` env. That matters because:

- `foli-decon`: R `4.2.3`, conda `r-music 0.2.0`
- old NEC notebook / reproduced sidecar: R `4.4.3`, GitHub `MuSiC 1.0.0`
  at `xuranw/MuSiC@f21fe67f5670d5e9fca0ad7550abaae3423eb59c`

Conda currently exposes `r-music 0.2.0`, while the exact NEC result requires the
newer GitHub MuSiC API:

```r
music_prop(bulk.mtx = ..., sc.sce = ...)
```

Trying to upgrade the main env to R 4.4 would conflict with the BayesPrism
constraint above. The compromise is a dedicated sidecar env:

```bash
mamba env create -f environment-r44.yml
mamba run -n foli-decon-r44 bash scripts/install_r44_sidecar.sh
```

The Python wrapper can call that sidecar with `r_env="foli-decon-r44"`.
For exact NEC replication it should also use a full SingleCellExperiment RDS
reference through `scrna_sce_rds_path`. Pre-subsetting the scRNA reference to the
Foli panel changes the MuSiC estimate even under MuSiC 1.0.0.

Current verified NEC hashes from the sidecar path:

- `Est.prop.weighted.tsv`: `468f0772bc3e45d51be29ce7a9370a0759e4be0e6a66fbff27cce1d5f1ca54c0`
- `Est.prop.allgene.tsv`: `e4be047da214ad4eaf7cc648fc3dd71834e3507751b9cf00f746df374ee05392`
- `Var.prop.tsv`: `9dcee2d83805f99dd7469250798a87b67190e9643f6f73e5fea1568f4d57de29`
- `Weight.gene.tsv`: `7d9be732e891d361094a2905afcf1997368394f48ac20e48e0ea54fa079b2153`
- `r.squared.full.tsv`: `0cfe88caa17e45ba79c582fb215e4a4a8f1b0ef6be0695d4e99054b3def917bc`
- `manifest.tsv`: `f662f9d5a98ccd45bcd59d369a03ab28b61f8e4c56584c090e676b3d70295d28`

Practical lesson: do not trade away BayesPrism by upgrading the main env. Use the
R 4.4 sidecar for MuSiC 1.0.0 and keep the main env R 4.2.

### xCell2 is the main non-conda straggler

`xCell2` is installed from GitHub, but its dependency tree is not allowed to
resolve itself from source. The current post-install call uses:

```r
remotes::install_github(
  "AlmogAngel/xCell2",
  upgrade = "never",
  dependencies = FALSE,
  build_vignettes = FALSE
)
```

This only works because `environment.yml` preinstalls the heavy imports through
conda, including packages such as `Rfast`, `ontologyIndex`, `minpack.lm`,
`singscore`, `BiocParallel`, and other Bioconductor infrastructure.

Current working snapshot:

- `xCell2 1.2.3`
- `Rfast 2.1.0`
- `xCell 1.1.0` loaded from R, with conda package `r-xcell`

Practical lesson: if `xCell2` is installed with `dependencies = TRUE`, R may try
to compile a large recursive dependency tree. That is exactly what the current
two-step install is designed to avoid.

### SCDC is conda-friendly, but R-version incompatible with the main env

SCDC is available as `r-scdc` from bioconda, but the available build requires
R `>=4.0,<4.1`. That conflicts with the main env's R 4.2 BayesPrism pin. The
wrapper uses the single-reference `SCDC::SCDC_prop` path.

Current repo change:

- `environment-scdc-r40.yml` defines `foli-decon-scdc-r40`
- wrapper: `run_scdc`
- native input class: `Biobase::ExpressionSet`

Observed dry-run:

- `mamba create -n foli-decon-scdc-probe -c conda-forge -c bioconda r-base=4.2 r-scdc --dry-run`
- result: unsatisfiable, because `r-scdc` requires R `<4.1`
- `mamba env create -f environment-scdc-r40.yml --dry-run`
- result: solved successfully as a sidecar with R `4.0.5`

Practical lesson: SCDC should be a sidecar method, but it can still be installed
with conda. Do not install it with `devtools::install_github("meichendong/SCDC")`
unless we are explicitly debugging the conda package.

### InstaPrism should stay in a sidecar

InstaPrism is GitHub-only and compiles Rcpp/RcppArmadillo code. Its dependency
list includes packages such as `caret`, `NMF`, `Rcpp`, and `RcppArmadillo`.
Those are all better installed by conda first.

Current repo compromise:

```bash
mamba env create -f environment-r44.yml
mamba run -n foli-decon-r44 bash scripts/install_r44_sidecar.sh
```

The install script checks conda-managed dependencies, then runs:

```r
remotes::install_github(
  "humengying0907/InstaPrism",
  dependencies = FALSE,
  upgrade = "never"
)
```

Practical lesson: this is the same source-install pattern as xCell2 and MuSiC
1.0.0, but InstaPrism is less central to the main env. Keep it isolated unless
we later find a conda build.

### CDSeq is conda-packaged, but only cleanly as an R 4.4 sidecar

CDSeq is available from bioconda as `r-cdseq`, which is exactly the kind of
package we prefer over source compilation. The catch is that the current
bioconda build targets R `>=4.4,<4.5` and pulls in a substantial stack including
Seurat, Harmony, Rcpp, RcppArmadillo, and Biobase. That conflicts with the main
R 4.2 stack chosen for BayesPrism.

Current repo compromise:

```bash
mamba env create -f environment-r44.yml
```

Use:

```python
run_cdseq(..., r_env="foli-decon-r44")
```

Practical lesson: CDSeq is not hard because it is unavailable; it is hard
because its clean conda build belongs to a newer R universe than the main env.
Keep it sidecar and avoid source-installing the full Seurat/Harmony tree into
`foli-decon`.

### Scaden and DISSECT can share one CPU TensorFlow sidecar

Scaden and DISSECT are both Python/TensorFlow deconvolution tools and are much
more compatible with each other than with the R-heavy main env. DISSECT is the
stricter package: its upstream requirements pin TensorFlow/Keras around 2.7 and
Python 3.9-era dependencies. Scaden only declares TensorFlow `>=2.0`, so it can
sit on the same TensorFlow 2.7 stack if installed without dependency upgrades.

Current repo compromise:

```bash
mamba env create -f environment-deep-cpu.yml
mamba run -n foli-decon-deep-cpu bash scripts/install_deep_cpu_tools.sh
```

The env intentionally keeps conda light: Python 3.9, Cython, compilers, and
`make`. The post-install script installs the old TensorFlow-era Python stack
with pip, then installs DISSECT and Scaden with `pip --no-deps` so pip does not
rewrite the chosen stack.

Observed solve behavior:

- first attempt: put TensorFlow/Scanpy/AnnData/h5py/pandas/scikit-learn pins
  directly in `environment-deep-cpu.yml`
- result: mamba did not reach a quick solution and was stopped after prolonged
  solving
- final attempt: keep only Python/build basics in conda and move the old Python
  stack to the post-install script
- result: `mamba env create -f environment-deep-cpu.yml --dry-run` solved as a
  46-package conda transaction

Practical lesson: this should remain a sidecar. The main package calls these
tools by subprocess, not import, so importing `foli_decon` does not load
TensorFlow and does not inherit old Python pins. CPU-only is simpler and more
portable; GPU support would need a separate CUDA/CuDNN-specific environment or
container.

### Requested expansion methods: compatibility sweep and smoke-tested sidecars

The following methods were explicitly requested for this branch and were checked
first by package/source probes and then, where feasible, by sidecar env creation
and small synthetic smoke tests.

The explicit compatibility matrix is maintained in
[`REQUESTED_TOOL_COMPATIBILITY_MATRIX.md`](REQUESTED_TOOL_COMPATIBILITY_MATRIX.md).

The latest sweep is recorded in
`scripts/check_requested_method_availability.sh` (2026-06-11):

- `r-tape`: `mamba search -c conda-forge -c bioconda r-tape` → no entries; real package route is PyPI `scTAPE`.
- `r-blade`: `mamba search -c conda-forge -c bioconda r-blade` → no entries; real package route is PyPI `BLADE-Deconvolution`.
- `r-mead`: `mamba search -c conda-forge -c bioconda r-mead` → no entries; the usable route is GitHub `DongyueXie/MEAD` plus GitHub dependency `mengyin/vashr`.
- `r-cellcode`: `mamba search -c conda-forge -c bioconda r-cellcode` → no entries; the usable route is GitHub `mchikina/CellCODE` with conda-managed `gplots` and `sva`.
- `r-dsa`: `mamba search -c conda-forge -c bioconda r-dsa` → no entries.
- `r-blue`: `mamba search -c conda-forge -c bioconda r-blue` → no entries; the real route is GitHub `SichenZhu/BLUE`, wrapped here as a CPU PyTorch/Scanpy sidecar pipeline adapter.
- `r-abis`: `mamba search -c conda-forge -c bioconda r-abis` → no entries; the usable route is the ABIS Shiny app source plus local signature resources.
- `r-estimate`: `mamba search -c conda-forge -c bioconda r-estimate` → no entries; the usable route is the original R-Forge `estimate 1.0.13` source package.
- `r-timer`: `mamba search -c conda-forge -c bioconda r-timer` → no entries.
- `r-timer3`: `mamba search -c conda-forge -c bioconda r-timer3` → no entries.
- `TIMER2`/`TIMER3`: blocked because the available route is web based, not a local package API.
- `bioconductor-psea`: `mamba search -c conda-forge -c bioconda bioconductor-psea` → available; `1.32.0` works with `r-base=4.2` and is now included in the main env spec.
- `r-consensustme`: `mamba search -c conda-forge -c bioconda r-consensustme` → available and integrated through the native `ConsensusTME` API.
- `r-mass` and `bioconductor-preprocesscore`: available and sufficient for the ABIS source logic; `preprocessCore` is only needed for ABIS microarray mode.
- `autogenes`: exists as `autogenes` Python package on bioconda and solves with
  `python=3.11`; feature selection and optional native deconvolution wrappers
  were added as a Python sidecar because the conda dependency set brings old
  AnnData-era packages.
- `markermap`: exists on PyPI as `markermap`; the package imports its comparison
  models at import time, including GitHub `iancovert/persist`, so the sidecar
  installs PyTorch/Lightning/Scanpy/COSG from conda and then installs
  `persist`, `lassonet`, and `markermap` with pip/source.
- `scGeneFit`: exists on PyPI as `scGeneFit`; install with `--no-deps` after
  conda supplies NumPy/SciPy/scikit-learn to avoid the deprecated PyPI
  `sklearn` dependency shim.
- `scTAPE`: exists on PyPI as `scTAPE==1.1.2`; wrapper added through a CPU
  PyTorch sidecar because naive pip resolution tries to pull CUDA Torch wheels.
- `BLADE-Deconvolution`: exists on PyPI as `BLADE-Deconvolution==0.0.7`;
  wrapper added through a small Python sidecar. The PyPI package named `BLADE`
  is unrelated iterator tooling and should not be used.
- `MEAD`: no conda package exists for MEAD or `vashr`; all compiled dependencies
  (`cluster`, `Rfast`, `ashr`, `SQUAREM`, `qvalue`, `truncnorm`, `mixsqp`,
  `quadprog`, `SingleCellExperiment`, and `SummarizedExperiment`) solve in
  `environment-r44.yml`. The post-install script installs only
  `mengyin/vashr` and `DongyueXie/MEAD` from GitHub with dependency upgrades
  disabled.
- `CellCODE`: no conda package exists, but `environment-r44.yml`
  solves with conda-managed `r-gplots` and `bioconductor-sva`; the post-install
  script installs only `mchikina/CellCODE` from GitHub.
- `BLUE`: source exists at `SichenZhu/BLUE`, and the wrapper deliberately treats
  the upstream numbered scripts as one pipeline-backed API. The adapter stages
  per-library h5ad files, a bulk TSV, cell-type mapping YAML, pipeline config,
  pseudobulk simulation, UNet training, and prediction into one `run_blue`
  call. The upstream `pyproject.toml` pins UV with CUDA 12.8 PyTorch wheels, so
  `environment-py311-cpu.yml` uses conda CPU PyTorch/Scanpy instead.
- PyPI name matches for `mead`, `blue`, `dsa`, and `timer` are not deconvolution
  packages. Observed dry-run metadata identified `mead` as Redis/aiohttp
  service tooling, `blue` as Python formatter tooling, CRAN `dsa` as daily time
  series seasonal adjustment, and PyPI `timer` as generic timing utilities.

Observed `mamba create --dry-run` checks:

- `mamba create --dry-run -n foli-decon-dry -c conda-forge -c bioconda r-base=4.2 r-consensustme=0.0.1.9000` → solved.
- `mamba create --dry-run -n foli-decon-dry -c conda-forge -c bioconda r-base=4.2 bioconductor-psea=1.32.0` → solved; wrapper added as PSEA marker reference-signal scoring.
- `available.packages(repos="http://r-forge.r-project.org")["estimate", ]` → `estimate 1.0.13`, `NeedsCompilation: no`; temp-library install loaded and exposed `estimateScore`, `filterCommonGenes`, and `outputGCT`.
- ABIS source probe downloaded `server.R`, `sigmatrixRNAseq.txt`, `sigmatrixMicro.txt`, and `target.txt` from `giannimonaco/ABIS`; the source formula uses `MASS::rlm` on intersected genes and multiplies coefficients by 100.
- `mamba create --dry-run -n foli-decon-autogenes-dry -c conda-forge -c bioconda python=3.11 autogenes=1.0.4` → solved.
- `mamba env create -f environment.yml --dry-run autogenes=1.0.4` → solved, but would pull conda `anndata 0.6.22.post1` before the pip `anndata>=0.10` install step, so the repo keeps AutoGeneS isolated in `environment-py311-cpu.yml`.
- `mamba env create --dry-run -f environment-r44.yml` → solved; merged R 4.4 sidecar for MuSiC 1.0.0/MuSiC2, MEAD, CellCODE, InstaPrism, and CDSeq.
- `mamba env create --dry-run -f environment-py311-cpu.yml` → solved; merged Python 3.11 CPU sidecar for AutoGeneS, BLADE, BLUE, TAPE, MarkerMap, and scGeneFit.
- `mamba env create --dry-run -f environment-scdc-r40.yml` → solved; SCDC remains separate because `r-scdc` requires R `<4.1`.
- `mamba env create --dry-run -f environment-deep-cpu.yml` → solved; Scaden and DISSECT remain together in the TensorFlow 2.7 sidecar.
- `scripts/check_requested_method_availability.sh` can be re-run to reproduce the exact same check with one command.

Observed real sidecar installs and smoke tests:

- `mamba env create -f environment-py311-cpu.yml` succeeded.
  `scripts/install_py311_cpu_sidecar.sh` installed `scTAPE==1.1.2`,
  `BLADE-Deconvolution==0.0.7`, `scGeneFit==1.0.2`, `lassonet==0.0.20`,
  GitHub `iancovert/persist`, and `markermap==1.0.3`, then cloned/validated
  upstream `SichenZhu/BLUE`.
- The merged Python sidecar import check loaded `autogenes`, `TAPE`,
  `BLADE_Deconvolution`, `scGeneFit`, `markermap`, `persist`, `torch 2.11.0`,
  and `scanpy 1.11.5`.
- Tiny synthetic wrapper smokes passed for `select_features_scgenefit`,
  `select_features_markermap`, `run_tape(..., epochs=1)`, `run_blade(...)`, and
  `run_blue(..., epochs=1)` using `python_env="foli-decon-py311-cpu"`.
  Outputs were written under `cache/merged_sidecar_smoke/`.
- TAPE required a runner-level PyTorch compatibility patch. Native TAPE saves
  and reloads a checkpoint created during the same temp run, but PyTorch 2.6+
  changed `torch.load` to `weights_only=True` by default. The sidecar runner now
  restores `weights_only=False` only inside `run_tape.py` before calling native
  TAPE. The wrapper also runs TAPE with its temp directory as `cwd` so native
  `model.pth` checkpoint side effects do not land in the repo root.
- MarkerMap smoke tests should avoid very small fixtures. Native MarkerMap
  asserts `max_epochs > min_epochs` and its BatchNorm layer can fail on a final
  one-sample training batch. The working smoke used 30 synthetic cells,
  `min_epochs=1`, `max_epochs=2`, and `batch_size=6`.
- `mamba env create -f environment-r44.yml` succeeded.
  `scripts/install_r44_sidecar.sh` installed GitHub MuSiC 1.0.0, `vashr`,
  MEAD, CellCODE, and InstaPrism over conda-managed R 4.4 dependencies. Load
  checks passed for `MuSiC 1.0.0`, `MEAD 1.0.2`, `CellCODE 0.99.0`,
  `InstaPrism 0.1.6`, and `CDSeq 1.0.9`.
- `mamba env create -f environment-deep-cpu.yml` succeeded.
  `scripts/install_deep_cpu_tools.sh` installed TensorFlow 2.7, DISSECT from
  GitHub, and `scaden==1.1.2` without compile failures. Import checks for
  TensorFlow, DISSECT, and `scaden --help` passed.
- Scaden wrapper smoke testing exposed two native I/O quirks: native pandas
  parsing tries to cast string cell IDs even with `index_col=0`, so the wrapper
  writes numeric temporary row indices; native `scaden simulate --prefix` writes
  merged h5ad paths relative to the current working directory, so the wrapper
  passes an absolute temp prefix.
- DISSECT wrapper smoke testing exposed an AnnData compatibility boundary:
  parent Python wrote newer nullable-string h5ad encodings that sidecar
  `anndata==0.8.0` could not read. The wrapper now stages TSV matrix/metadata
  and the sidecar runner writes h5ad using DISSECT's own AnnData version.

Practical lesson for this expansion pass:
- Keep requested unsupported methods out of the unified env until a clear, stable
  conda or source strategy exists. PSEA, ABIS, ESTIMATE, AutoGeneS, TAPE,
  BLADE, BLUE, MEAD, CellCODE, MarkerMap, scGeneFit, Scaden, and DISSECT are the
  exceptions from this sweep because they have compatible package routes,
  narrow source/resource wrappers, or smoke-tested sidecar environments.
- For methods not in channels, explicit failure guidance is now surfaced from
  `run_deconvolution` so users get actionable diagnostics before any runtime
  work starts.

### MCP-counter hides a runtime network dependency

`r-mcpcounter` itself installs cleanly through conda, but the native
`MCPcounter::MCPcounter.estimate` defaults fetch marker tables from GitHub:

- `Signatures/probesets.txt`
- `Signatures/genes.txt`

That means a "local" MCP-counter run can fail in a network-restricted runtime
even when the package is installed. The wrapper now exposes `probesets_path` and
`genes_path` so marker tables can be supplied from local files.

Observed smoke tests:

- `run_quantiseq` completed with installed native `quantiseqr` and a
  `(2, 11)` output matrix whose first sample summed to `1.0`
- `run_mcp_counter` failed without local marker tables because
  `raw.githubusercontent.com` could not resolve
- `run_mcp_counter` succeeded with a local four-row `probesets.tsv`, producing
  a `(3, 2)` score matrix

Practical lesson: MCP-counter is not an install blocker, but it is a
reproducibility blocker unless marker tables are cached or passed explicitly.

### AnnotationHub remains a small Bioconductor post-install

`AnnotationHub` is installed by `BiocManager` in `scripts/install_r_packages.R`.
It is not part of the heavy deconvolution solver stack, but `xCell2` expects it.
Keeping this as a post-install step has been simpler than forcing it into the
conda solve for the R 4.2 stack.

Current working snapshot:

- `AnnotationHub 3.6.0`

Practical lesson: this is acceptable as a small post-install, but it should not
become a precedent for installing the full tool stack from R.

### CIBERSORTx is not an R dependency problem

CIBERSORTx is isolated in the `cibersortx/fractions` container. The main hurdle
was rootless Podman behavior, not package dependency resolution.

On this machine, plain `podman run` tried to use `pasta`, which was unavailable,
and `crun` also failed for this image. The working wrapper default is:

```bash
podman run --runtime=runc --network=slirp4netns ...
```

Current working snapshot:

- `podman 5.5.0`
- `runc 1.4.1`
- `slirp4netns 1.1.8`

Practical lesson: keep CIBERSORTx containerized and keep the Podman runtime args
explicit in the wrapper. This avoids mixing CIBERSORTx concerns into the R env.

### dtangle installs cleanly, but is a runtime/data hurdle

`dtangle` was not a major environment creation blocker. The working env installs
it directly from conda:

- conda package: `r-dtangle 2.0.9`

The hurdle showed up during panel-data benchmarking instead. On the NEC panel
benchmarks, every `dtangle` row failed after the wrapper detected non-finite
predictions:

- `benchmark_out_20260322_panel_run3/combined_metrics.tsv`
- `benchmark_out_20260322_keysplit_full/combined_metrics.tsv`

This is different from the BayesPrism/xCell2 dependency issues. It means the
package can be installed, loaded, and called, but the current default wrapper
settings are not robust for these Foli-seq panel simulations.

Likely causes to investigate later:

- too few informative marker genes after panel intersection
- marker selection producing degenerate reference sets
- transformed panel matrices containing rows with little or no usable variation
- `n_markers=50` being too aggressive for small panel/reference combinations

Practical lesson: keep `dtangle` in the unified env, but treat it as a method
that needs panel-specific tuning before using its benchmark metrics. Good next
experiments are lower `n_markers`, stricter gene filtering before `dtangle`, and
explicit checks for all-NA/non-finite estimates from the native result.

## What Not To Do

Avoid these install patterns unless intentionally testing a new environment:

- Do not install the full R stack with `install.packages()`,
  `BiocManager::install()`, or `remotes::install_github(..., dependencies = TRUE)`.
- Do not casually update R beyond 4.2 without checking `r-bayesprism`.
- Do not let `xCell2` upgrade dependencies during GitHub install.
- Do not install GitHub MuSiC 1.0.0 into the main R 4.2 env for exact NEC
  replication; use the `foli-decon-r44` sidecar.
- Do not treat CIBERSORTx as something to install into R or Python.

## Current Compromise

The compromise is:

- one unified conda env for Python and almost all R deconvolution tools
- R pinned to 4.2 for BayesPrism compatibility
- one merged R 4.4 sidecar for MuSiC 1.0.0/MuSiC2, MEAD, CellCODE, InstaPrism, and CDSeq
- an R 4.0 sidecar for SCDC because the conda build requires R `<4.1`
- a Python 3.9/TensorFlow 2.7 CPU sidecar for Scaden and DISSECT
- one merged Python 3.11 CPU sidecar for AutoGeneS, BLADE, BLUE, TAPE, MarkerMap, and scGeneFit
- conda-managed compiled dependencies wherever possible
- a tiny post-install for `AnnotationHub` and `xCell2`
- local MCP-counter marker tables for offline/reproducible score runs
- CIBERSORTx isolated in a container, called through the same Python API

This is still simpler than per-tool environments for every method, while avoiding
the most painful recursive source compilation path.

## Relearning Checklist

If this has to be rebuilt from scratch, test in this order:

1. `mamba env create -f environment.yml`
2. `mamba run -n foli-decon Rscript scripts/install_r_packages.R`
3. `mamba env create -f environment-r44.yml`
4. `mamba run -n foli-decon-r44 bash scripts/install_r44_sidecar.sh`
5. `mamba env create -f environment-scdc-r40.yml`
6. `mamba env create -f environment-deep-cpu.yml`
7. `mamba run -n foli-decon-deep-cpu bash scripts/install_deep_cpu_tools.sh`
8. `mamba env create -f environment-py311-cpu.yml`
9. `mamba run -n foli-decon-py311-cpu bash scripts/install_py311_cpu_sidecar.sh`
10. `mamba run -n foli-decon pytest -q`
11. A direct `run_cibersortx(...)` smoke test
12. A `run_mcp_counter(...)` smoke test with local marker tables
13. Small `select_features_scgenefit` and `select_features_markermap` smoke tests with tiny synthetic scRNA data
14. Small `run_cdseq`, `run_tape`, `run_blade`, `run_scaden`, and `run_dissect` smoke tests with tiny iteration counts
15. A small benchmark with `music,cibersortx,xcell2`

If step 1 fails, inspect the R 4.2 conda solve first. If step 2 fails, inspect
`xCell2` and `AnnotationHub` before changing unrelated packages.
