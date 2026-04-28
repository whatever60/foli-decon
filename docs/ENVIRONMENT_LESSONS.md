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
mamba env create -f environment-music-r44.yml
mamba run -n foli-decon-music-r44 Rscript scripts/install_music_r44.R
```

The Python wrapper can call that sidecar with `r_env="foli-decon-music-r44"`.
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
  replication; use the `foli-decon-music-r44` sidecar.
- Do not treat CIBERSORTx as something to install into R or Python.

## Current Compromise

The compromise is:

- one unified conda env for Python and almost all R deconvolution tools
- R pinned to 4.2 for BayesPrism compatibility
- an R 4.4 sidecar only for MuSiC 1.0.0 exact-replication workflows
- conda-managed compiled dependencies wherever possible
- a tiny post-install for `AnnotationHub` and `xCell2`
- CIBERSORTx isolated in a container, called through the same Python API

This is still simpler than per-tool environments for every method, while avoiding
the most painful recursive source compilation path.

## Relearning Checklist

If this has to be rebuilt from scratch, test in this order:

1. `mamba env create -f environment.yml`
2. `mamba run -n foli-decon Rscript scripts/install_r_packages.R`
3. `mamba env create -f environment-music-r44.yml`
4. `mamba run -n foli-decon-music-r44 Rscript scripts/install_music_r44.R`
5. `mamba run -n foli-decon pytest -q`
6. A direct `run_cibersortx(...)` smoke test
7. A small benchmark with `music,cibersortx,xcell2`

If step 1 fails, inspect the R 4.2 conda solve first. If step 2 fails, inspect
`xCell2` and `AnnotationHub` before changing unrelated packages.
