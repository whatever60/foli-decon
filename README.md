# foli-decon

`foli-decon` is a Python package that wraps multiple cell type deconvolution tools behind one interface for Foli-seq style panel RNA-seq data.

Foli-seq characteristics assumed by this repo:
- amplicon panel expression (~500 genes)
- UMI count matrix as mixture input
- frequent gene mismatch between mixture and reference

This repo standardizes three things:
1. input normalization/transform steps (for example `counts -> CPM -> log1p`)
2. gene intersection between mixture and reference before tool execution
3. output shape (`samples x cell_types`)

## Supported tools

| Tool | Wrapper | Expected reference | Typical transform | Extra required inputs |
|---|---|---|---|---|
| xCell | `run_xcell` | none (built-in signatures) | `cpm` | none |
| xCell2 | `run_xcell2` | trained `xCell2` object (`.rds`) | `cpm` | `xcell2_object_path` |
| EPIC | `run_epic` | built-in or custom signature (`genes x cell_types`) | `cpm` or `tpm` | optional `tumor` |
| dtangle | `run_dtangle` | signature (`genes x cell_types`) | `cpm_log1p` | optional `n_markers` |
| BayesPrism | `run_bayesprism` | scRNA counts (`genes x cells`) + cell types | `raw` | optional cell states, tumor key, outlier controls |
| DeconRNASeq | `run_deconrnaseq` | signature (`genes x cell_types`) | `cpm` | optional `checksig` |
| CIBERSORTx | `run_cibersortx` | signature (`genes x cell_types`) | `raw` | `username`, `token`, container runtime |
| Bisque | `run_bisque` | scRNA counts + cell types + donor IDs | `raw` | `batch_ids` |
| DWLS | `run_dwls` | signature (`genes x cell_types`) | `cpm` | optional submethod |
| MuSiC | `run_music` | scRNA counts + cell types + donor IDs, or full SCE RDS | `raw` | `batch_ids`, or `scrna_sce_rds_path` + metadata columns |

## Install

### 1. Create the conda/mamba environment

```bash
mamba env create -f environment.yml
mamba activate foli-decon
```

Note: `r-base` is pinned to `4.2` in this unified stack because the conda `r-bayesprism` builds currently target R `<4.3`.
This environment also installs `podman` for `CIBERSORTx` container execution.

### 2. Install non-conda stragglers (currently xCell2)

```bash
Rscript scripts/install_r_packages.R
```

This script intentionally does not install the full tool stack from source.
It expects core R packages (including heavy `xCell2` imports like `Rfast`) to already be present from `environment.yml`.
It installs `AnnotationHub` from Bioconductor (not conda-available for this R pin) and then installs only `xCell2` from GitHub (`dependencies = FALSE`).

### 3. Optional MuSiC R 4.4 sidecar for exact MuSiC 1.0.0 runs

The unified `foli-decon` env stays on R 4.2 for `BayesPrism`. For MuSiC workflows that need the newer `MuSiC 1.0.0` SingleCellExperiment API or bit-identical replication of the NEC notebooks, create the sidecar env:

```bash
mamba env create -f environment-music-r44.yml
mamba run -n foli-decon-music-r44 Rscript scripts/install_music_r44.R
```

This pins `xuranw/MuSiC` to commit `f21fe67f5670d5e9fca0ad7550abaae3423eb59c` under R 4.4.3 without upgrading the main env.

Environment build lessons and known dependency compromises are recorded in
[`docs/ENVIRONMENT_LESSONS.md`](docs/ENVIRONMENT_LESSONS.md).

## Quick start

```python
import pandas as pd
from foli_decon import run_deconvolution

mixture = pd.read_csv("foli_counts.tsv", sep="\t", index_col=0)
signature = pd.read_csv("reference_signature.tsv", sep="\t", index_col=0)

result = run_deconvolution(
    tool="deconrnaseq",
    mixture=mixture,
    signature=signature,
    mixture_transform="cpm",
    signature_transform="cpm",
)

print(result.proportions.head())
```

## Two-split benchmark workflow

The package includes a simulation benchmark flow to select deconvolution tools for a specific scRNA-seq reference.

Workflow:
1. split scRNA-seq cells into `reference` and `evaluation` partitions by cell type
2. build reference objects from the reference split
3. simulate pseudobulk mixtures from evaluation cells with known truth
4. run each deconvolution tool
5. score predictions against truth

Implemented simulation scenarios:
- `dirichlet_uniform`
- `one_type_dominant`
- `one_cell_type_dominate`

Implemented benchmark APIs:
- `run_two_split_benchmark_from_adata(...)`
- `run_two_split_benchmark_from_h5ad(...)`

Main outputs:
- `metrics` table with one row per `dataset x tool x scenario`
- primary metrics: `mae`, `rmse`, `js_divergence`, `cosine_similarity`
- additional divergence metrics: `kl_divergence_truth_to_pred`, `kl_divergence_pred_to_truth`, `kl_divergence_symmetric`, `hellinger_distance`, `total_variation_distance`, `bray_curtis_distance`
- overlap metadata: `n_shared_samples`, `n_shared_cell_types`
- status fields: `status`, `error_message`

CLI helper script:

```bash
mamba run -n foli-decon python scripts/run_two_split_benchmark.py \
  --adata-paths \
    /home/ubuntu/dev/20251201_nec/dropbox/20250911_NEC_SI_atlas/smallintestine_nec_processed.h5ad \
    /home/ubuntu/dev/20251201_nec/dropbox/20250911_NEC_SI_atlas/smallintestine_nec_processed_subset_foli_deconvolution.h5ad \
  --cell-type-col celltype \
  --batch-col donor \
  --count-layer raw_counts \
  --output-dir /home/ubuntu/dev/foli_decon/benchmark_out
```

Notes:
- `xcell2` needs `--xcell2-object-path`.
- use `--xcell2-min-shared-genes` for panel inputs when overlap with the xCell2 reference is low.
- `cibersortx` needs `--cibersortx-username` and `--cibersortx-token`.
- panel benchmarking is supported through `--panel-genes-path`.
- benchmark scoring is skipped for `xcell2` because it returns non-compositional scores.
- default split is key-based (`key_stratified`) using `batch_col`; override with `--split-strategy random_cell_type` if needed.
- set `--split-key-col donor` (or sample/batch key) to control which key defines reference/eval partitioning.

## Input conventions

All wrappers use these conventions:
- matrices are `pandas.DataFrame`
- genes in rows
- samples or cells in columns
- `cell_types` and `batch_ids` are vectors aligned to scRNA matrix columns
- wrappers write temporary TSV files and call backend engines through subprocesses

## API reference by tool

Unified dispatcher:

```python
run_deconvolution(tool: str, **kwargs) -> DeconvolutionResult
```

Per-tool wrappers:

```python
run_xcell(
    mixture: pd.DataFrame,
    transform: str = "cpm",
    arrays: bool = False,
    expected_cell_types: list[str] | None = None,
) -> DeconvolutionResult

run_xcell2(
    mixture: pd.DataFrame,
    xcell2_object_path: str,
    transform: str = "cpm",
    min_shared_genes: float = 0.2,
    raw_scores: bool = False,
    spillover: bool = True,
    spillover_alpha: float = 0.5,
) -> DeconvolutionResult

run_epic(
    mixture: pd.DataFrame,
    signature: pd.DataFrame | None = None,
    scrna_counts: pd.DataFrame | None = None,
    cell_types: pd.Series | None = None,
    mixture_transform: str = "cpm",
    signature_transform: str = "cpm",
    gene_lengths: pd.Series | None = None,
    tumor: bool = False,
    scale_mrna: bool = True,
) -> DeconvolutionResult

run_dtangle(
    mixture: pd.DataFrame,
    signature: pd.DataFrame | None = None,
    scrna_counts: pd.DataFrame | None = None,
    cell_types: pd.Series | None = None,
    mixture_transform: str = "cpm_log1p",
    signature_transform: str = "cpm_log1p",
    gene_lengths: pd.Series | None = None,
    n_markers: int = 50,
    marker_method: str = "ratio",
) -> DeconvolutionResult

run_deconrnaseq(
    mixture: pd.DataFrame,
    signature: pd.DataFrame | None = None,
    scrna_counts: pd.DataFrame | None = None,
    cell_types: pd.Series | None = None,
    mixture_transform: str = "cpm",
    signature_transform: str = "cpm",
    gene_lengths: pd.Series | None = None,
    checksig: bool = False,
    use_scale: bool = True,
) -> DeconvolutionResult

run_dwls(
    mixture: pd.DataFrame,
    signature: pd.DataFrame | None = None,
    scrna_counts: pd.DataFrame | None = None,
    cell_types: pd.Series | None = None,
    mixture_transform: str = "cpm",
    signature_transform: str = "cpm",
    gene_lengths: pd.Series | None = None,
    dwls_submethod: str = "DampenedWLS",
) -> DeconvolutionResult

run_music(
    mixture: pd.DataFrame,
    scrna_counts: pd.DataFrame,
    cell_types: pd.Series,
    batch_ids: pd.Series,
    mixture_transform: str = "raw",
    scrna_transform: str = "raw",
    select_ct: list[str] | None = None,
) -> DeconvolutionResult

run_bisque(
    mixture: pd.DataFrame,
    scrna_counts: pd.DataFrame,
    cell_types: pd.Series,
    batch_ids: pd.Series,
    mixture_transform: str = "raw",
    scrna_transform: str = "raw",
    use_overlap: bool = False,
    old_cpm: bool = True,
) -> DeconvolutionResult

run_bayesprism(
    mixture: pd.DataFrame,
    scrna_counts: pd.DataFrame,
    cell_types: pd.Series,
    mixture_transform: str = "raw",
    scrna_transform: str = "raw",
    cell_states: pd.Series | None = None,
    tum_key: str | None = None,
    update_gibbs: bool = True,
    n_cores: int = 1,
    which_theta: str = "final",
    state_or_type: str = "type",
    outlier_cut: float = 1.0,
    outlier_fraction: float = 1.0,
    pseudo_min: float = 1e-8,
) -> DeconvolutionResult

run_cibersortx(
    mixture: pd.DataFrame,
    username: str,
    token: str,
    signature: pd.DataFrame | None = None,
    scrna_counts: pd.DataFrame | None = None,
    cell_types: pd.Series | None = None,
    mixture_transform: str = "raw",
    signature_transform: str = "raw",
    container_runtime: str = "podman",
    image: str = "cibersortx/fractions",
    label: str | None = None,
    container_runtime_args: list[str] | None = None,
) -> DeconvolutionResult
```

### xCell (`run_xcell`)

- Engine: R via subprocess (`Rscript`), calling `immunedeconv::deconvolute_xcell`.
- Required args: `mixture` (`pd.DataFrame`, genes x samples).
- Optional args: `transform`, `arrays`, `expected_cell_types`.
- Special features: optional array mode and expected cell type subset.
- Native backend input expectation: in-memory expression matrix (R matrix); no signature file required.
- Wrapper preprocessing: validate numeric matrix, normalize gene names (`.` -> `-`), apply transform.

### xCell2 (`run_xcell2`)

- Engine: R via subprocess (`Rscript`), calling `xCell2Analysis` from package `xCell2` with helper script `xCell2Analysis.R`.
- Required args: `mixture` (`pd.DataFrame`), `xcell2_object_path` (`str`, path to trained `.rds` xCell2 object).
- Optional args: `transform`, `min_shared_genes`, `raw_scores`, `spillover`, `spillover_alpha`.
- Special features: supports lineage/context-specific trained xCell2 reference objects.
- Native backend input expectation: expression matrix + serialized xCell2 object loaded with `readRDS`.
- Wrapper preprocessing: validate matrix, normalize genes, apply transform.

### EPIC (`run_epic`)

- Engine: R via subprocess, calling `immunedeconv::deconvolute_epic` or `immunedeconv::deconvolute_epic_custom`.
- Required args: `mixture`.
- Optional args: custom reference via `signature` or (`scrna_counts` + `cell_types`), transforms, `gene_lengths` (for TPM modes), `tumor`, `scale_mrna`.
- Special features: built-in EPIC reference mode or custom signature mode.
- Native backend input expectation: expression matrix; custom mode needs `signature_matrix` and `signature_genes`.
- Wrapper preprocessing: validate, normalize genes, optional signature construction from scRNA by cell-type averaging, transform, mixture/signature gene intersection.

### dtangle (`run_dtangle`)

- Engine: R via subprocess, calling `dtangle::dtangle`.
- Required args: `mixture` plus `signature` or (`scrna_counts` + `cell_types`).
- Optional args: transforms, `gene_lengths`, `n_markers`, `marker_method`.
- Special features: marker selection controls (`n_markers`, `marker_method`).
- Native backend input expectation: RNA-seq-like matrices with samples x genes and references x genes on a compatible scale.
- Wrapper preprocessing: optional signature from scRNA, normalize genes, transform (default `cpm_log1p`), intersection, transpose to dtangle orientation.

### DeconRNASeq (`run_deconrnaseq`)

- Engine: R via subprocess, solved with `limSolve::lsei` in the R wrapper (not direct call to the `DeconRNASeq` package API).
- Required args: `mixture` plus `signature` or (`scrna_counts` + `cell_types`).
- Optional args: transforms, `gene_lengths`, `checksig`, `use_scale`.
- Special features: constrained solution with non-negativity + sum-to-one enforced in `lsei`; `checksig` currently accepted but ignored.
- Native backend input expectation: bulk sample vector and signature matrix in a linear system.
- Wrapper preprocessing: optional signature build from scRNA, normalize genes, transform (default `cpm`), intersection.

### DWLS (`run_dwls`)

- Engine: R via subprocess, calling `DWLS::solveDampenedWLS` / `DWLS::solveSVR` / `DWLS::solveOLS`.
- Required args: `mixture` plus `signature` or (`scrna_counts` + `cell_types`).
- Optional args: transforms, `gene_lengths`, `dwls_submethod` in `{"OLS","SVR","DampenedWLS"}`.
- Special features: selectable internal solver backend.
- Native backend input expectation: signature matrix and per-sample bulk vector with matched genes.
- Wrapper preprocessing: optional signature build from scRNA, normalize genes, transform (default `cpm`), intersection.

### MuSiC (`run_music`)

- Engine: R via subprocess, calling `MuSiC::music_prop`.
- Required args: `mixture` plus either (`scrna_counts`, `cell_types`, `batch_ids`) or `scrna_sce_rds_path`.
- Optional args: `mixture_transform`, `scrna_transform`, `select_ct`, `scrna_sce_rds_path`, `cell_type_column`, `batch_id_column`, `r_env`, `r_output_path`, `extra_output_dir`.
- Special features: explicit donor/sample ID (`batch_ids`) support to model cross-subject variability.
- Special features: can call the R 4.4 MuSiC sidecar via `r_env="foli-decon-music-r44"` and a full SingleCellExperiment RDS reference via `scrna_sce_rds_path`.
- Native backend input expectation: main matrix path uses `Biobase::ExpressionSet` for bulk and scRNA with `cellType` and `sampleID` metadata; SCE/RDS path uses `MuSiC 1.0.0` native `bulk.mtx` + `sc.sce`.
- Wrapper preprocessing: matrix path validates metadata lengths, normalizes genes, intersects mixture/scRNA, transforms, and writes aligned vectors; SCE/RDS path transforms only the mixture and lets MuSiC handle common genes inside the full reference.

### Bisque (`run_bisque`)

- Engine: R via subprocess, calling `BisqueRNA::ReferenceBasedDecomposition`.
- Required args: `mixture`, `scrna_counts`, `cell_types`, `batch_ids`.
- Optional args: `mixture_transform`, `scrna_transform`, `use_overlap`, `old_cpm`.
- Special features: donor-aware reference decomposition and overlap toggles.
- Native backend input expectation: `ExpressionSet` objects with scRNA metadata fields `cellType` and `batchId`.
- Wrapper preprocessing: validate/align metadata, normalize genes, intersection, transform, serialize TSV inputs.

### BayesPrism (`run_bayesprism`)

- Engine: R via subprocess, calling `BayesPrism::new.prism`, `BayesPrism::run.prism`, and `BayesPrism::get.fraction`.
- Required args: `mixture`, `scrna_counts`, `cell_types`.
- Optional args: `mixture_transform`, `scrna_transform`, `cell_states`, `tum_key`, `update_gibbs`, `n_cores`, `which_theta`, `state_or_type`, `outlier_cut`, `outlier_fraction`, `pseudo_min`.
- Special features: supports cell state labels, tumor key designation, and outlier controls useful for panel data.
- Native backend input expectation: count matrices (`input.type = "count.matrix"`), with reference as cells x genes and mixture as samples x genes.
- Wrapper preprocessing: normalize genes, intersection, transform (defaults `raw`/`raw` to preserve count semantics), write `cell_types` and optional `cell_states`.

### CIBERSORTx (`run_cibersortx`)

- Engine: container CLI via subprocess (`podman`/`docker`) calling `cibersortx/fractions`.
- Required args: `mixture`, `username`, `token`, plus one reference source: `signature` or (`scrna_counts` + `cell_types`).
- Optional args: `mixture_transform`, `signature_transform`, `container_runtime`, `image`, `label`, `container_runtime_args`.
- Special features: local container compute with CIBERSORTx credentialed execution; supports scRNA-to-signature conversion.
- Native backend input expectation: mounted TSV files (`mixture.tsv`, `sigmatrix.tsv`) and CLI flags; output is `CIBERSORTx_Results.txt` (or label-specific filename).
- Wrapper preprocessing: normalize genes, transform, intersection, write TSVs, then post-process output by dropping non-proportion columns (`Correlation`, `RMSE`, `P-value`).
- Podman default runtime args are `--runtime=runc --network=slirp4netns` to avoid rootless `pasta` network failures.

## Output data type and format

Common return type for every wrapper:

- `DeconvolutionResult.tool`: `str`
- `DeconvolutionResult.proportions`: `pd.DataFrame` with shape `samples x outputs` and numeric values (typically float)
- `DeconvolutionResult.metadata`: `dict[str, str | int | float | bool]` containing wrapper/runtime settings

Output semantics by tool:

- `xCell`: enrichment scores per sample and cell type. These are not compositional fractions and are not constrained to sum to 1.
- `xCell2`: enrichment-like scores per sample and cell type. This is score output, not ratio output; behavior depends on `raw_scores` and spillover options.
- `EPIC`: fraction-like output matrix from EPIC/immunedeconv. Columns can include non-immune/other categories depending on reference mode.
- `dtangle`: estimated mixture fractions from `result$estimates`.
- `DeconRNASeq` wrapper: constrained coefficients from `limSolve::lsei` (non-negative, sum-to-one constraints are enforced in this wrapper).
- `DWLS`: estimated composition matrix from selected solver (`OLS`, `SVR`, or `DampenedWLS`).
- `MuSiC`: weighted fraction estimate (`Est.prop.weighted`) is returned. When `extra_output_dir` is supplied for the SCE/RDS path, supplementary native outputs (`Est.prop.allgene`, `Weight.gene`, `Var.prop`, `r.squared.full`, and a manifest) are also written to disk.
- `Bisque`: only `bulk.props` (deconvolved proportions) is returned.
- `BayesPrism`: returns only `theta` fractions (`get.fraction` output when `update_gibbs=True`, otherwise initial posterior theta). Output can be type-level or state-level depending on `state_or_type`.
- `CIBERSORTx`: returns only cell-type proportion columns. Supplemental columns from native output (`Correlation`, `RMSE`, `P-value`) are removed.

Current limitation:

- Wrappers intentionally normalize to one main matrix output. Most backend-specific diagnostics, uncertainty measures, and intermediate objects are not yet exposed in `DeconvolutionResult`.

## Gene mismatch behavior by method

Important correction: methods in this package do **not** use gene union for deconvolution.  
Most paths use gene intersection, and some tools apply additional internal overlap checks.

General wrapper behavior:
- for explicit matrix references (`signature`, `scrna_counts`), wrappers usually normalize names then intersect genes before backend calls
- `xCell` and `xCell2` are special cases because references are internal (`xCell`) or serialized model objects (`xCell2`)

Method-by-method mismatch handling:

| Method | Wrapper-level handling | Backend-level handling | Panel-data implication |
|---|---|---|---|
| xCell | no explicit wrapper reference intersection | `xCell::rawEnrichmentAnalysis` intersects with xCell gene universe and fails when shared genes `< 5000` | high risk for panel data (~500 genes); usually not tolerant |
| xCell2 | wrapper passes `min_shared_genes` (default `0.2`) to `xCell2Analysis` | `xCell2Analysis` intersects with `getGenesUsed(xcell2object)`; stops if shared fraction `< minSharedGenes`, warns if `< 0.85` | high risk if reference is whole-transcriptome and mixture is panel; prefer panel-matched xCell2 object |
| EPIC | custom-reference path intersects `mixture/signature` before call | `EPIC::EPIC` intersects bulk/reference genes internally; warns when common genes are low (`<2000`); stops if signature genes in common are fewer than reference cell types | moderate risk with very small overlap or too many cell types |
| dtangle | wrapper intersects `mixture/signature` before call | marker selection runs on provided genes; no hard overlap threshold in `dtangle::dtangle` itself | moderate risk; low overlap means weaker marker set |
| DeconRNASeq wrapper (`limSolve::lsei`) | wrapper intersects `mixture/signature` before call | constrained least-squares solve; no explicit overlap threshold | moderate risk when genes are too few relative to number of cell types |
| DWLS | wrapper intersects `mixture/signature`; R layer intersects again | solver runs on intersected genes; no explicit high-overlap threshold in wrapper | moderate risk; marker-poor panels can destabilize estimates |
| MuSiC | matrix path intersects `mixture/scrna` before call; SCE/RDS path keeps the full reference and lets MuSiC intersect internally | `MuSiC::music_prop` checks common genes and can stop with `Too few common genes!`; full-reference SCE/RDS mode can produce different estimates than pre-subset panel references | usually tolerant, but for exact replication prefer full-reference SCE/RDS mode with the MuSiC R 4.4 sidecar |
| Bisque | wrapper intersects `mixture/scrna` before call | `BisqueRNA` overlap utilities require non-zero overlap and can stop if filtering leaves zero genes | moderate risk if filtering removes most panel genes |
| BayesPrism | wrapper intersects `mixture/scrna` before call | `BayesPrism::new.prism` intersects again; stops at zero shared genes; warns when shared genes `< 100` | moderate risk for very tiny panels; ~500-gene panels generally satisfy this warning threshold |
| CIBERSORTx | wrapper intersects `mixture/signature` before container run | native fractions image overlap minimum is not explicitly surfaced by wrapper; output computed on intersected genes | unknown-to-moderate risk; empirical wrapper runs succeeded with 20 and 120 shared genes in this repo |

Why this matters for Foli-seq:
- panel assays often keep only a small subset of transcriptome genes
- methods with hard overlap requirements (`xCell`, strict `xCell2` settings) are most likely to fail
- for other methods, runs may complete but precision can degrade as overlap/marker informativeness drops

## Backends

- R-backed wrappers call `src/foli_decon/r_scripts/run_tools.R` via `Rscript` subprocess
- CIBERSORTx wrapper calls containerized `cibersortx/fractions` via `podman`/`docker` subprocess
- current implementation does not use `rpy2`; R communication is file + subprocess based
- shared preprocess logic lives in `src/foli_decon/preprocess.py`

## Notes on tool-specific extras

- `xCell2` can require lineage-specific references (for example immuno-oncology); pass the correct `.rds` object.
- `MuSiC` and `Bisque` require donor/sample IDs (`batch_ids`) for reference cells.
- `BayesPrism` can use subtype labels (`cell_states`) in addition to primary `cell_types`.
- `BayesPrism` exposes `outlier_cut`, `outlier_fraction`, and `pseudo_min` for panel-size tuning.
- `CIBERSORTx` requires web credentials and container runtime access.

## CIBERSORTx wrapper usage

`run_cibersortx` is available through both direct import and `run_deconvolution(tool="cibersortx", ...)`.

```python
import pandas as pd
from foli_decon import run_cibersortx

mixture = pd.read_csv("foli_counts.tsv", sep="\t", index_col=0)
signature = pd.read_csv("sigmatrix.tsv", sep="\t", index_col=0)

result = run_cibersortx(
    mixture=mixture,
    signature=signature,
    username="your_cibersortx_account",
    token="your_cibersortx_token",
    container_runtime="podman",
)

print(result.proportions.head())
```

Notes:
- CIBERSORTx container runs require `username` and `token` (password is not used by the fractions image).
- The wrapper aligns mixture/signature genes before writing container input files.

## Current strategy for dependency conflicts

Primary path:
1. single environment `foli-decon`
2. conda/mamba packages first for all available decon tools
3. small post-install script only for packages missing from conda (currently `xCell2`)

Fallback path if a method remains incompatible in the unified env:
1. keep Python wrapper API stable
2. run that method in a dedicated container
3. return the same `samples x cell_types` output format
