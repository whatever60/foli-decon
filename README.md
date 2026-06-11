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
| quanTIseq | `run_quantiseq` | none (built-in immune signatures) | `cpm` | optional `tumor`, `arrays`, `scale_mrna` |
| MCP-counter | `run_mcp_counter` | none (built-in/native marker tables) | `cpm` | optional local `probesets_path`/`genes_path` to avoid native GitHub downloads |
| PSEA | `run_psea` | marker sets keyed by cell type | `cpm` | `marker_sets`; returns marker reference signals, not fractions |
| CellCODE | `run_cellcode` | pure cell-type signature or scRNA-derived pseudobulk signature | `cpm_log1p` | bulk sample `groups`; returns surrogate proportion variable scores, not fractions |
| ABIS | `run_abis` | ABIS RNA-seq/microarray signature matrices | RNA-seq default `cpm`; native expects TPM-like RNA-seq input | run `scripts/download_abis_resources.py` or pass `signature_path`; output is native percent rescaled to fractions by default |
| ESTIMATE | `run_estimate` | none (built-in stromal/immune signatures from R-Forge package) | `cpm` | install source package via `scripts/install_r_packages.R`; returns scores, not cell-type fractions |
| CDSeq | `run_cdseq` | optional signature/scRNA reference for component labeling; can run reference-free | `raw` | `cell_type_number` unless a reference is supplied; sidecar env recommended |
| BLADE | `run_blade` | signature mean/sd or scRNA counts + cell types | mixture `raw`, reference `cpm -> log1p` | Python sidecar recommended; signature SD is required or estimated |
| BLUE | `run_blue` | scRNA counts + cell types, plus optional donor/sample IDs | `raw` | CPU PyTorch sidecar plus upstream BLUE repo; returns proportions and optional ctGEP h5ad |
| xCell2 | `run_xcell2`, `train_xcell2_reference` | trained `xCell2` object (`.rds`) or scRNA counts + cell types | `cpm` for mixture, `raw` for scRNA training | `xcell2_object_path`, or `scrna_counts` + `cell_types` |
| EPIC | `run_epic` | built-in or custom signature (`genes x cell_types`) | `cpm` or `tpm` | optional `tumor` |
| dtangle | `run_dtangle` | signature (`genes x cell_types`) | `cpm_log1p` | optional `n_markers` |
| BayesPrism | `run_bayesprism` | scRNA counts (`genes x cells`) + cell types | `raw` | optional cell states, tumor key, outlier controls |
| DeconRNASeq | `run_deconrnaseq` | signature (`genes x cell_types`) | `cpm` | optional `checksig` |
| CIBERSORTx | `run_cibersortx` | signature (`genes x cell_types`) | `raw` | `username`, `token`, container runtime |
| Bisque | `run_bisque` | scRNA counts + cell types + donor IDs | `raw` | `batch_ids` |
| DWLS | `run_dwls` | signature (`genes x cell_types`) | `cpm` | optional submethod |
| MuSiC | `run_music` | scRNA counts + cell types + donor IDs, or full SCE RDS | `raw` | `batch_ids`, or `scrna_sce_rds_path` + metadata columns |
| MuSiC2 | `run_music2` | scRNA counts + cell types + donor IDs, or full SCE RDS | `raw` | separate `control_mixture` and `case_mixture` |
| MEAD | `run_mead` | scRNA counts + cell types + donor IDs | `raw` | sidecar env recommended; optional bulk sample `groups`; returns proportions plus standard-error uncertainty |
| SCDC | `run_scdc` | scRNA counts + cell types + donor IDs | `raw` | `batch_ids`; optional `ct_sub` |
| InstaPrism | `run_instaprism` | scRNA counts + cell types; optional cell states | `raw` | optional `cell_states`; sidecar env recommended |
| TAPE | `run_tape` | scRNA counts + cell types | `raw` | CPU PyTorch sidecar; trains native TAPE autoencoder |
| Scaden | `run_scaden` | scRNA counts + cell types, or an existing trained Scaden model | `raw` | CPU TensorFlow sidecar; trains native Scaden ensemble |
| DISSECT | `run_dissect` | scRNA counts + cell types; optional batch IDs | mixture `cpm`, scRNA `raw` | CPU TensorFlow sidecar; semi-supervised training |
| AutoGeneS | `select_features_autogenes`, `run_autogenes` | signature or scRNA-derived pseudobulk signature | `cpm` | Python sidecar recommended; feature selection first, optional native deconvolution workflow |
| MarkerMap | `select_features_markermap` | scRNA counts + cell types | `cpm_log1p` | CPU PyTorch/Lightning feature-selector sidecar; returns selected genes, not fractions |
| scGeneFit | `select_features_scgenefit` | scRNA counts + cell types | `cpm_log1p` | lightweight Python feature selector; returns selected genes, not fractions |

## Requested additions reviewed

The following methods were evaluated for availability during setup. Supported sidecar methods in this table have also had at least a small synthetic smoke test when noted in `docs/ENVIRONMENT_LESSONS.md`.
The explicit compatibility matrix is maintained in [`docs/REQUESTED_TOOL_COMPATIBILITY_MATRIX.md`](docs/REQUESTED_TOOL_COMPATIBILITY_MATRIX.md).

| Method | Package probe | Dry-run outcome | Status in repo |
|---|---|---|---|
| TAPE | `r-tape`; real package `scTAPE` | R alias absent; PyPI `scTAPE==1.1.2` and `environment-py311-cpu.yml` solve | supported through CPU PyTorch sidecar (`run_tape`) |
| BLADE | `r-blade`; real package `BLADE-Deconvolution` | R alias absent; PyPI `BLADE-Deconvolution==0.0.7` and `environment-py311-cpu.yml` solve | supported through Python sidecar (`run_blade`) |
| MEAD | `r-mead`; GitHub `DongyueXie/MEAD`; GitHub dependency `mengyin/vashr`; PyPI `mead` is unrelated | no conda package for MEAD/vashr, but all compiled dependencies solve in `environment-r44.yml` | supported through R sidecar plus GitHub post-install (`run_mead`) |
| ABIS | `r-abis` absent; GitHub Shiny app source/resources available | conda alias absent; local source logic works with `MASS::rlm` and ABIS signature resources | supported through Shiny-source-derived local wrapper (`run_abis`) |
| ESTIMATE | `r-estimate` absent; R-Forge `estimate 1.0.13` available | conda alias absent; no-compile source package installs from R-Forge | supported through original `estimate::estimateScore` (`run_estimate`) |
| CONSENSUS_TME | `r-consensustme` | available and integrated (`run_consensus_tme`) | supported through native `ConsensusTME` |
| TIMER | `r-timer` | no conda package found | blocked; local native route unresolved |
| TIMER2/TIMER3 | web tools | web based only | blocked for local wrappers |
| BLUE | `r-blue`; PyPI `blue`; GitHub `SichenZhu/BLUE` | conda/PyPI names are unrelated; real code is a UV/PyTorch pipeline with YAML configs and numbered scripts; CPU conda sidecar solves | supported through CPU PyTorch sidecar pipeline adapter (`run_blue`) |
| PSEA | `bioconductor-psea` | available and integrated (`run_psea`) | supported as marker reference-signal scores |
| DSA | `r-dsa`; CRAN/PyPI `dsa` | no conda deconvolution package found; CRAN/PyPI names are unrelated | blocked until a stable local package/source route is chosen |
| CellCODE | `r-cellcode` / `bioconductor-cellcode`; GitHub `mchikina/CellCODE` | no conda package found; conda deps solve in `environment-r44.yml`; GitHub package returns SPV scores, not fractions | supported as score output through R sidecar (`run_cellcode`) |
| AutoGeneS | `autogenes` | solves with `python=3.11`; full main-env solve also works but pulls old conda AnnData before pip upgrades | supported through Python sidecar (`select_features_autogenes`, `run_autogenes`) |
| MarkerMap | PyPI `markermap`; GitHub `Computational-Morphogenomics-Group/MarkerMap` | conda base stack solves in `environment-py311-cpu.yml`; `markermap` and GitHub `iancovert/persist` install as pip/source layer | supported as feature selector (`select_features_markermap`) |
| scGeneFit | PyPI `scGeneFit`; GitHub `solevillar/scGeneFit-python` | conda base stack solves in `environment-py311-cpu.yml`; pip package installs without compiling | supported as feature selector (`select_features_scgenefit`) |

Reproduce the full method availability sweep with:

```bash
scripts/check_requested_method_availability.sh
```

## Install

### 1. Create the conda/mamba environment

```bash
mamba env create -f environment.yml
mamba activate foli-decon
```

Note: `r-base` is pinned to `4.2` in this unified stack because the conda `r-bayesprism` builds currently target R `<4.3`.
This environment also installs `podman` for `CIBERSORTx` container execution.

### 2. Install non-conda stragglers (`estimate`, `xCell2`)

```bash
Rscript scripts/install_r_packages.R
```

This script intentionally does not install the full tool stack from source.
It expects core R packages (including heavy `xCell2` imports like `Rfast`) to already be present from `environment.yml`.
It installs `AnnotationHub` from Bioconductor (not conda-available for this R pin), installs the no-compile R-Forge `estimate` package, and then installs only `xCell2` from GitHub (`dependencies = FALSE`).

### 3. Download optional ABIS source resources

ABIS is a Shiny app, not an R package. The local wrapper uses the app's published signature matrices and `MASS::rlm` logic.

```bash
python scripts/download_abis_resources.py
```

By default `run_abis` looks for:

- `resources/abis/rnaseq/sigmatrixRNAseq.txt`
- `resources/abis/microarray/sigmatrixMicro.txt`
- `resources/abis/microarray/target.txt`

### 4. Optional merged R 4.4 sidecar

The unified `foli-decon` env stays on R 4.2 for `BayesPrism`. Tools that need R 4.4 are compressed into one sidecar:

- MuSiC2 and exact MuSiC 1.0.0 workflows
- MEAD
- CellCODE
- InstaPrism
- CDSeq

```bash
mamba env create -f environment-r44.yml
mamba run -n foli-decon-r44 bash scripts/install_r44_sidecar.sh
```

Use these tools with `r_env="foli-decon-r44"`.

### 5. Optional SCDC R 4.0 sidecar

SCDC is also isolated because the available conda `r-scdc` build requires R 4.0, which conflicts with the main R 4.2/BayesPrism env:

```bash
mamba env create -f environment-scdc-r40.yml
```

Use `run_scdc(..., r_env="foli-decon-scdc-r40")` after creating this env.

### 6. Optional CPU deep-learning sidecar for Scaden and DISSECT

Scaden and DISSECT can share a CPU-only TensorFlow 2.7 sidecar. The conda env only provides Python/build basics; the post-install script installs the old TensorFlow-era Python wheels and then installs Scaden/DISSECT without dependency upgrades:

```bash
mamba env create -f environment-deep-cpu.yml
mamba run -n foli-decon-deep-cpu bash scripts/install_deep_cpu_tools.sh
```

Use `run_scaden(..., python_env="foli-decon-deep-cpu")` or `run_dissect(..., python_env="foli-decon-deep-cpu")`.

### 7. Optional merged Python 3.11 CPU sidecar

Most Python sidecar tools can share one CPU PyTorch/Scanpy env:

- AutoGeneS
- BLADE
- BLUE
- TAPE
- MarkerMap
- scGeneFit

```bash
mamba env create -f environment-py311-cpu.yml
mamba run -n foli-decon-py311-cpu bash scripts/install_py311_cpu_sidecar.sh
```

Use these tools with `python_env="foli-decon-py311-cpu"`.
By default `scripts/install_py311_cpu_sidecar.sh` clones upstream BLUE into `resources/blue/BLUE`. Set `FOLI_DECON_BLUE_REPO=/path/to/BLUE` or pass `blue_repo_path=...` to use a different checkout.

Environment build lessons and known dependency compromises are recorded in
[`docs/ENVIRONMENT_LESSONS.md`](docs/ENVIRONMENT_LESSONS.md).
The main env and active sidecar YAMLs were dry-run solved on 2026-06-11.

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

print(result.proportion.head())
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
- `xcell2` can use `--xcell2-object-path`; custom xCell2 references can also be trained first with `train_reference(tool="xcell2", ...)`.
- use `--xcell2-min-shared-genes` for panel inputs when overlap with the xCell2 reference is low.
- `cibersortx` needs `--cibersortx-username` and `--cibersortx-token`.
- panel benchmarking is supported through `--panel-genes-path`.
- proportion-error metrics are skipped for `xcell2` because it returns non-compositional scores; score-association metrics are reported when labels overlap.
- default split is key-based (`key_stratified`) using `batch_col`; override with `--split-strategy random_cell_type` if needed.
- set `--split-key-col donor` (or sample/batch key) to control which key defines reference/eval partitioning.
- `autogenes`, `blade`, `blue`, `cdseq`, `cibersortx`, `dissect`, `instaprism`, `music2`, `psea`, `scaden`, `scdc`, and `tape` are not included in the default CLI tool list because they need credentials, marker sets, sidecar envs, heavier training, or non-single-mixture inputs; pass them explicitly when configured.

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
train_reference(tool: str, **kwargs) -> ReferenceTrainingResult
select_features(tool: str, **kwargs) -> FeatureSelectionResult
```

Per-tool wrappers:

```python
run_xcell(
    mixture: pd.DataFrame,
    transform: str = "cpm",
    arrays: bool = False,
    expected_cell_types: list[str] | None = None,
) -> DeconvolutionResult

run_quantiseq(
    mixture: pd.DataFrame,
    transform: str = "cpm",
    tumor: bool = False,
    arrays: bool = False,
    scale_mrna: bool = True,
) -> DeconvolutionResult

run_mcp_counter(
    mixture: pd.DataFrame,
    transform: str = "cpm",
    feature_types: str = "HUGO_symbols",
    probesets_path: str | None = None,
    genes_path: str | None = None,
) -> DeconvolutionResult

run_psea(
    mixture: pd.DataFrame,
    marker_sets: dict[str, list[str] | list[list[str]]],
    transform: str = "cpm",
    sample_subset: list[str] | None = None,
    target_mean: float = 1.0,
) -> DeconvolutionResult

run_cellcode(
    mixture: pd.DataFrame,
    signature: pd.DataFrame | None = None,
    scrna_counts: pd.DataFrame | None = None,
    cell_types: pd.Series | None = None,
    groups: pd.Series | None = None,
    mixture_transform: str = "cpm_log1p",
    signature_transform: str = "cpm_log1p",
    tag_cutoff: float = 2.0,
    max_markers: int | None = None,
    method: str = "mixed",
    mix_par: float = 0.3,
    r_env: str | None = None,
) -> DeconvolutionResult

run_autogenes(
    mixture: pd.DataFrame,
    signature: pd.DataFrame | None = None,
    scrna_counts: pd.DataFrame | None = None,
    cell_types: pd.Series | None = None,
    mixture_transform: str = "cpm",
    signature_transform: str = "cpm",
    gene_lengths: pd.Series | None = None,
    python_env: str | None = None,
    ngen: int = 50,
    nfeatures: int = 500,
    mode: str = "fixed",
    seed: int = 0,
    selection_index: int = 0,
    model: str = "nusvr",
) -> DeconvolutionResult

select_features_autogenes(
    signature: pd.DataFrame | None = None,
    scrna_counts: pd.DataFrame | None = None,
    cell_types: pd.Series | None = None,
    signature_transform: str = "cpm",
    gene_lengths: pd.Series | None = None,
    python_env: str | None = None,
    ngen: int = 50,
    nfeatures: int = 500,
    mode: str = "fixed",
    seed: int = 0,
    selection_index: int = 0,
) -> FeatureSelectionResult

select_features_markermap(
    scrna_counts: pd.DataFrame,
    cell_types: pd.Series,
    scrna_transform: str = "cpm_log1p",
    gene_lengths: pd.Series | None = None,
    python_env: str | None = None,
    nfeatures: int = 50,
    hidden_layer_size: int | None = None,
    z_size: int = 16,
    batch_size: int = 64,
    loss_tradeoff: float = 0.0,
    train_fraction: float = 0.8,
    min_epochs: int = 2,
    max_epochs: int = 10,
    auto_lr: bool = False,
    seed: int = 0,
) -> FeatureSelectionResult

select_features_scgenefit(
    scrna_counts: pd.DataFrame,
    cell_types: pd.Series,
    scrna_transform: str = "cpm_log1p",
    gene_lengths: pd.Series | None = None,
    python_env: str | None = None,
    nfeatures: int = 50,
    method: str = "centers",
    epsilon: float = 1.0,
    sampling_rate: float = 1.0,
    n_neighbors: int = 3,
    max_constraints: int = 1000,
    redundancy: float = 0.01,
    seed: int = 0,
) -> FeatureSelectionResult

run_blade(
    mixture: pd.DataFrame,
    signature: pd.DataFrame | None = None,
    signature_sd: pd.DataFrame | None = None,
    scrna_counts: pd.DataFrame | None = None,
    cell_types: pd.Series | None = None,
    mixture_transform: str = "raw",
    signature_transform: str = "cpm",
    signature_is_log: bool = False,
    python_env: str | None = None,
    n_rep: int = 3,
    n_jobs: int = 1,
    n_rep_final: int = 3,
) -> DeconvolutionResult

run_blue(
    mixture: pd.DataFrame,
    scrna_counts: pd.DataFrame,
    cell_types: pd.Series,
    batch_ids: pd.Series | None = None,
    mixture_transform: str = "raw",
    scrna_transform: str = "raw",
    python_env: str | None = None,
    blue_repo_path: str | None = None,
    celltype_mapping: dict[str, list[str]] | None = None,
    input_genes: list[str] | None = None,
    output_genes: list[str] | None = None,
    n_training_samples: int = 3000,
    n_validation_samples: int = 300,
    cells_per_sample: int = 500,
    epochs: int = 50,
    batch_size: int = 64,
    extra_output_dir: str | None = None,
) -> DeconvolutionResult

run_cdseq(
    mixture: pd.DataFrame,
    signature: pd.DataFrame | None = None,
    scrna_counts: pd.DataFrame | None = None,
    cell_types: pd.Series | None = None,
    mixture_transform: str = "raw",
    signature_transform: str = "raw",
    gene_lengths: pd.Series | None = None,
    cell_type_number: int | list[int] | None = None,
    beta: float | list[float] | None = 0.5,
    alpha: float = 5.0,
    mcmc_iterations: int = 700,
    dilution_factor: float = 1.0,
    gene_subset_size: int | None = None,
    block_number: int = 1,
    cpu_number: int | None = 1,
    r_env: str | None = None,
) -> DeconvolutionResult

run_xcell2(
    mixture: pd.DataFrame,
    xcell2_object_path: str | None = None,
    scrna_counts: pd.DataFrame | None = None,
    cell_types: pd.Series | None = None,
    xcell2_reference_output_path: str | None = None,
    transform: str = "cpm",
    scrna_transform: str = "raw",
    return_signatures: bool = False,
    min_shared_genes: float = 0.2,
    raw_scores: bool = False,
    spillover: bool = True,
    spillover_alpha: float = 0.5,
) -> DeconvolutionResult

train_xcell2_reference(
    scrna_counts: pd.DataFrame,
    cell_types: pd.Series,
    output_path: str,
    return_signatures: bool = False,
) -> ReferenceTrainingResult

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

run_music2(
    control_mixture: pd.DataFrame,
    case_mixture: pd.DataFrame,
    scrna_counts: pd.DataFrame | None = None,
    cell_types: pd.Series | None = None,
    batch_ids: pd.Series | None = None,
    control_transform: str = "raw",
    case_transform: str = "raw",
    scrna_transform: str = "raw",
    select_ct: list[str] | None = None,
    scrna_sce_rds_path: str | None = None,
    r_env: str | None = None,
) -> DeconvolutionResult

run_mead(
    mixture: pd.DataFrame,
    scrna_counts: pd.DataFrame,
    cell_types: pd.Series,
    batch_ids: pd.Series,
    mixture_transform: str = "raw",
    scrna_transform: str = "raw",
    select_ct: list[str] | None = None,
    marker_genes: list[str] | None = None,
    filter_gene: bool = False,
    groups: pd.Series | None = None,
    calc_var: bool = True,
    use_qp: bool = False,
    r_env: str | None = None,
) -> DeconvolutionResult

run_scdc(
    mixture: pd.DataFrame,
    scrna_counts: pd.DataFrame,
    cell_types: pd.Series,
    batch_ids: pd.Series,
    mixture_transform: str = "raw",
    scrna_transform: str = "raw",
    ct_sub: list[str] | None = None,
    r_env: str | None = None,
) -> DeconvolutionResult

run_instaprism(
    mixture: pd.DataFrame,
    scrna_counts: pd.DataFrame,
    cell_types: pd.Series,
    mixture_transform: str = "raw",
    scrna_transform: str = "raw",
    cell_states: pd.Series | None = None,
    r_env: str | None = None,
) -> DeconvolutionResult

run_tape(
    mixture: pd.DataFrame,
    scrna_counts: pd.DataFrame,
    cell_types: pd.Series,
    mixture_transform: str = "raw",
    scrna_transform: str = "raw",
    python_env: str | None = None,
    variance_threshold: float = 0.98,
    scaler: str = "mms",
    datatype: str = "counts",
    epochs: int = 128,
) -> DeconvolutionResult

run_scaden(
    mixture: pd.DataFrame,
    scrna_counts: pd.DataFrame | None = None,
    cell_types: pd.Series | None = None,
    mixture_transform: str = "raw",
    scrna_transform: str = "raw",
    python_env: str | None = None,
    model_dir: str | None = None,
    retrain: bool = True,
    n_training_samples: int = 1000,
    cells_per_sample: int = 100,
    train_steps: int = 5000,
) -> DeconvolutionResult

run_dissect(
    mixture: pd.DataFrame,
    scrna_counts: pd.DataFrame,
    cell_types: pd.Series,
    batch_ids: pd.Series | None = None,
    mixture_transform: str = "cpm",
    scrna_transform: str = "raw",
    python_env: str | None = None,
    n_training_samples: int = 1000,
    cells_per_sample: int = 500,
    train_steps: int = 5000,
    n_models: int = 5,
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

- Engine: R via subprocess (`Rscript`), calling `xCell::xCellAnalysis`.
- Required args: `mixture` (`pd.DataFrame`, genes x samples).
- Optional args: `transform`, `arrays`, `expected_cell_types`.
- Special features: optional array mode and expected cell type subset.
- Native backend input expectation: in-memory expression matrix (R matrix); no signature file required.
- Wrapper preprocessing: validate numeric matrix, normalize gene names (`.` -> `-`), apply transform.

### quanTIseq (`run_quantiseq`)

- Engine: R via subprocess, calling `quantiseqr::run_quantiseq`.
- Required args: `mixture` (`pd.DataFrame`, genes x samples).
- Optional args: `transform`, `tumor`, `arrays`, `scale_mrna`.
- Special features: tumor-mode marker filtering and optional correction for cell-type-specific mRNA content.
- Native backend input expectation: non-log expression matrix with gene symbols as row names.
- Wrapper preprocessing: validate, normalize gene names, apply transform (default `cpm`).
- Output: fraction-like matrix with immune cell types plus an `other` column; rows should sum to 1 in normal runs.

### MCP-counter (`run_mcp_counter`)

- Engine: R via subprocess, calling `MCPcounter::MCPcounter.estimate` directly.
- Required args: `mixture` (`pd.DataFrame`, genes/probes x samples).
- Optional args: `transform`, `feature_types`, `probesets_path`, `genes_path`.
- Special features: returns abundance scores, not proportions; supports `HUGO_symbols`, `ENTREZ_ID`, and `affy133P2_probesets`.
- Native backend input expectation: expression matrix plus marker tables. If marker tables are not supplied, native MCP-counter defaults try to download them from GitHub.
- Wrapper preprocessing: validate, normalize gene names, apply transform (default `cpm`), and pass local marker tables when provided.
- Output: score matrix; values are not constrained to sum to 1.

### PSEA (`run_psea`)

- Engine: R via subprocess, calling `PSEA::marker`.
- Required args: `mixture` and `marker_sets` (`dict[str, list[str] | list[list[str]]]`) where keys are cell type names.
- Optional args: `transform`, `sample_subset`, `target_mean`.
- Special features: marker groups can represent multiple probes/features for one transcript by passing nested lists.
- Native backend input expectation: expression matrix with marker feature IDs present in row names; `sampleSubset` is integer column indices, which the wrapper can derive from supplied sample names.
- Wrapper preprocessing: validate, normalize gene names, apply transform (default `cpm`), normalize marker gene names, and write a marker-set TSV.
- Output: score matrix of PSEA population reference signals. These are relative marker signals, not compositional fractions.

### CellCODE (`run_cellcode`)

- Engine: R via subprocess, calling `CellCODE::tagData` and `CellCODE::getAllSPVs`.
- Required args: `mixture`, bulk sample `groups`, plus `signature` or (`scrna_counts` + `cell_types`).
- Optional args: transforms, `gene_lengths`, `tag_cutoff`, `max_markers`, `ref_mean`, `method`, `mix_par`, `r_env`, `extra_output_dir`.
- Special features: estimates surrogate proportion variables while modeling group effects; this is useful for latent cell-composition signals, not direct fraction estimation.
- Native backend input expectation: expression matrix genes x samples, group labels with at least two levels, and a marker tag matrix genes x cell types built from pure cell-type expression.
- Wrapper preprocessing: optional scRNA-to-signature aggregation, normalize genes, transform (default `cpm_log1p`), intersect mixture/signature genes, generate CellCODE marker tags in R.
- Output: sample x cell type SPV score matrix in `score`. Optional extras write `data_tag.tsv`, `spv.tsv`, and a manifest.

### ABIS (`run_abis`)

- Engine: R via subprocess, using the ABIS Shiny app source logic (`MASS::rlm`) rather than `immunedeconv`.
- Required args: `mixture` (`pd.DataFrame`, genes x samples). Also provide `signature_path` or run `scripts/download_abis_resources.py` so the default local ABIS signature files exist.
- Optional args: `transform`, `technology` (`rnaseq` or `microarray`), `signature`, `signature_path`, `target_path`, `gene_lengths`, `rescale_percent`, `r_env`, `output_path`, `log_path`.
- Special features: supports RNA-seq and microarray ABIS signatures; native fit is unconstrained robust linear regression, so negative coefficients are possible.
- Native backend input expectation: ABIS RNA-seq expects TPM-like matrix values with gene symbols; microarray mode additionally quantile-normalizes to the ABIS target distribution.
- Wrapper preprocessing: validate, normalize gene names, apply transform (default `cpm`), pass mixture/signature to R, and let ABIS source logic intersect genes.
- Output: native ABIS percent estimates are returned as `proportion` after division by 100 by default. Set `rescale_percent=False` to keep native percent scale.

### ESTIMATE (`run_estimate`)

- Engine: R via subprocess, calling the original R-Forge `estimate::estimateScore`.
- Required args: `mixture` (`pd.DataFrame`, genes x samples).
- Optional args: `transform`, `platform`, `gene_lengths`, `r_env`, `output_path`, `log_path`.
- Special features: tumor purity is returned only for native platforms that support it, especially `affymetrix`; RNA-seq-style runs usually return stromal, immune, and ESTIMATE scores.
- Native backend input expectation: GCT v1.2 file with gene symbols as `NAME` rows and samples as columns.
- Wrapper preprocessing: validate, normalize gene names, apply transform (default `cpm`), write native GCT input, parse native GCT output into sample x score TSV.
- Output: `score` matrix with `StromalScore`, `ImmuneScore`, `ESTIMATEScore`, and sometimes `TumorPurity`; not cell-type proportions.

### AutoGeneS (`select_features_autogenes`, `run_autogenes`)

- Engine: Python sidecar subprocess, calling `autogenes.init`, `autogenes.optimize`, `autogenes.select`, and optionally `autogenes.deconvolve`.
- Required args for feature selection: `signature` or (`scrna_counts` + `cell_types`).
- Required args for native deconvolution workflow: `mixture` plus `signature` or (`scrna_counts` + `cell_types`).
- Optional args: transforms, `gene_lengths`, `python_env`, optimizer controls (`ngen`, `nfeatures`, `mode`, `seed`, `selection_index`), deconvolution `model`, coefficient clipping/normalization, and output paths.
- Special features: primarily performs multi-objective marker/feature selection; `run_autogenes` additionally exposes native AutoGeneS deconvolution with models `nusvr`, `nnls`, and `linear`.
- Native backend input expectation: reference profile as cell types x genes and bulk matrix as samples x genes.
- Wrapper preprocessing: optional scRNA-to-signature aggregation, gene-name normalization, transform, mixture/signature intersection, and TSV/JSON payload creation for the sidecar runner.
- Output: `select_features_autogenes` returns `FeatureSelectionResult.selected_genes`; `run_autogenes` returns a sample x cell type coefficient matrix plus `DeconvolutionResult.selected_features`. By default the deconvolution wrapper clips negative coefficients and row-normalizes them into proportions.

### MarkerMap (`select_features_markermap`)

- Engine: Python sidecar subprocess, calling native `markermap.vae_models.MarkerMap` and `train_model`.
- Required args: `scrna_counts` (`pd.DataFrame`, genes x cells) and `cell_types` (`pd.Series`, aligned to scRNA columns).
- Optional args: `scrna_transform`, `gene_lengths`, `python_env`, `nfeatures`, VAE/training controls (`hidden_layer_size`, `z_size`, `batch_size`, `loss_tradeoff`, `train_fraction`, `min_epochs`, `max_epochs`, `auto_lr`, `seed`), and output paths.
- Special features: supervised nonlinear marker selection; `loss_tradeoff=0.0` makes the default run supervised-only for cell-type discrimination.
- Native backend input expectation: AnnData cells x genes with `.obs["cell_type"]`; the sidecar builds this object from staged TSV files.
- Wrapper preprocessing: validates scRNA counts and labels, normalizes gene names, applies transform, writes genes x cells TSV plus label TSV, and lets the sidecar build native AnnData.
- Output: `FeatureSelectionResult.selected_genes` plus `ranking` with `rank`, `gene`, and native feature index. No proportions or scores are returned.

### scGeneFit (`select_features_scgenefit`)

- Engine: Python sidecar subprocess, calling native `scGeneFit.functions.get_markers`.
- Required args: `scrna_counts` (`pd.DataFrame`, genes x cells) and `cell_types` (`pd.Series`, aligned to scRNA columns).
- Optional args: `scrna_transform`, `gene_lengths`, `python_env`, `nfeatures`, native LP controls (`method`, `epsilon`, `sampling_rate`, `n_neighbors`, `max_constraints`, `redundancy`), `seed`, and output paths.
- Special features: label-aware linear-programming marker selection; `method="centers"` is the fastest/stablest native default.
- Native backend input expectation: dense cells x genes NumPy matrix plus one label per cell.
- Wrapper preprocessing: validates scRNA counts and labels, normalizes gene names, applies transform, converts labels to stable integer category codes in the sidecar, and writes TSV/JSON payloads.
- Output: `FeatureSelectionResult.selected_genes` plus `ranking` with `rank`, `gene`, and native feature index. No proportions or scores are returned.

### BLADE (`run_blade`)

- Engine: Python sidecar subprocess, importing `BLADE-Deconvolution`.
- Required args: `mixture` plus either (`signature`, optional `signature_sd`) or (`scrna_counts` + `cell_types`).
- Optional args: transforms, `signature_is_log`, `signature_sd_floor`, `python_env`, BLADE hyperparameter grids, repeat counts, `n_jobs`, and `extra_output_dir`.
- Special features: jointly estimates composition and cell-type-specific expression while modeling reference variability.
- Native backend input expectation: signature mean matrix `genes x cell_types` in log scale, signature SD matrix `genes x cell_types`, and linear-scale bulk matrix `genes x samples`.
- Wrapper preprocessing: normalizes genes, transforms and intersects mixture/reference genes, builds log-scale mean/SD from scRNA when requested, and clips SD values to avoid zero-variance failures.
- Output: `final_obj.ExpF(final_obj.Beta)` as sample x cell type proportions. Optional extras include selected hyperparameters and run IDs.

### CDSeq (`run_cdseq`)

- Engine: R via subprocess, calling `CDSeq::CDSeq`.
- Required args: `mixture`; additionally `cell_type_number` when no reference is supplied.
- Optional args: `signature`, or (`scrna_counts` + `cell_types`) for reference-assisted component labeling; `gene_lengths`, `beta`, `alpha`, `mcmc_iterations`, `dilution_factor`, `gene_subset_size`, `block_number`, `cpu_number`, `r_env`, `extra_output_dir`.
- Special features: reference-free complete deconvolution; can estimate latent cell type count when `cell_type_number` is a vector; returns cell-type-specific GEPs as optional extras.
- Native backend input expectation: raw RNA-seq read count matrix with genes x samples; optional reference GEP with the same genes.
- Wrapper preprocessing: validates, normalizes gene names, optionally constructs a pseudobulk reference from scRNA by summing cells per type, intersects mixture/reference genes, and writes TSV inputs.
- Output: sample x component matrix from `t(estProp)`. Without a reference, columns are latent `CDSeq_estimated_cell_type_*`; with a reference, CDSeq tries to assign/reference-label components.

### xCell2 (`run_xcell2`)

- Engine: R via subprocess (`Rscript`), calling `xCell2Analysis` from package `xCell2` with helper script `xCell2Analysis.R`.
- Training engine: R via subprocess, calling `xCell2::xCell2Train`.
- Required args for scoring a pre-trained object: `mixture` (`pd.DataFrame`) and `xcell2_object_path` (`str`, path to trained `.rds` xCell2 object).
- Required args for end-to-end train-then-score: `mixture`, `scrna_counts`, and `cell_types`; pass `xcell2_reference_output_path` to persist the trained object.
- Separate training API: `train_reference(tool="xcell2", ...)` or direct `train_xcell2_reference(...)`.
- Optional args: `transform`, `scrna_transform`, `min_shared_genes`, `raw_scores`, `spillover`, `spillover_alpha`, `return_signatures`, `min_pb_cells`, `min_pb_samples`, `min_sc_genes`, `xcell2_workers`.
- Special features: supports lineage/context-specific trained xCell2 reference objects and signature-only references via `return_signatures=True`.
- Native backend input expectation: expression matrix + serialized xCell2 object loaded with `readRDS`; training uses scRNA-seq count matrix + labels.
- Wrapper preprocessing: validate matrices, normalize genes, apply requested transforms.

### EPIC (`run_epic`)

- Engine: R via subprocess, calling `EPIC::EPIC`.
- Required args: `mixture`.
- Optional args: custom reference via `signature` or (`scrna_counts` + `cell_types`), transforms, `gene_lengths` (for TPM modes), `tumor`, `scale_mrna`.
- Special features: built-in EPIC reference mode or custom signature mode.
- Native backend input expectation: expression matrix; custom mode passes an EPIC reference list with `refProfiles` and `sigGenes`.
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
- Special features: can call the R 4.4 MuSiC sidecar via `r_env="foli-decon-r44"` and a full SingleCellExperiment RDS reference via `scrna_sce_rds_path`.
- Native backend input expectation: main matrix path uses `Biobase::ExpressionSet` for bulk and scRNA with `cellType` and `sampleID` metadata; SCE/RDS path uses `MuSiC 1.0.0` native `bulk.mtx` + `sc.sce`.
- Wrapper preprocessing: matrix path validates metadata lengths, normalizes genes, intersects mixture/scRNA, transforms, and writes aligned vectors; SCE/RDS path transforms only the mixture and lets MuSiC handle common genes inside the full reference.

### MuSiC2 (`run_music2`)

- Engine: R via subprocess, calling `MuSiC::music2_prop_t_statistics`.
- Required args: `control_mixture`, `case_mixture`, plus either (`scrna_counts`, `cell_types`, `batch_ids`) or `scrna_sce_rds_path`.
- Optional args: transforms, `select_ct`, SCE metadata columns, `r_env`, `extra_output_dir`, and MuSiC2 iteration/DE-gene cutoffs.
- Special features: condition-aware deconvolution; estimates case samples while detecting/removing cell-type-specific DE genes between control and case bulk groups.
- Native backend input expectation: two bulk matrices, a `SingleCellExperiment`, cluster column, donor/sample column, and selected cell types.
- Wrapper preprocessing: matrix path intersects control/case/scRNA genes; SCE/RDS path intersects control/case genes and leaves scRNA intersection to MuSiC2.
- Output: `Est.prop` proportions. When `extra_output_dir` is supplied, convergence, iteration count, removed DE genes, and non-converged sample IDs are written when present.

### MEAD (`run_mead`)

- Engine: R via subprocess, calling `MEAD::MEAD` in the MEAD sidecar env.
- Required args: `mixture`, `scrna_counts`, `cell_types`, `batch_ids`.
- Optional args: transforms, `select_ct`, `marker_genes`, `gene_thresh`, quantile filters, `filter_gene`, `hc_type`, `centering_xy`, `nfold`, bulk sample `groups`, `calc_var`, `use_qp`, `r_env`, `extra_output_dir`.
- Special features: donor-aware scRNA reference through MEAD's `individual` metadata; optional bulk sample group labels can produce group-difference diagnostics; `uncertainty` stores standard errors when `calc_var=True`.
- Native backend input expectation: bulk count matrix `genes x samples` and a `SingleCellExperiment` reference with `colData(ref)$cell_type` and `colData(ref)$individual`.
- Wrapper preprocessing: validates metadata, normalizes genes, applies transforms, intersects mixture/scRNA genes, and builds the `SingleCellExperiment` inside R.
- Output: `p_hat` proportions as sample x cell type. Optional extras can include `p_hat`, `p_hat_se`, confidence intervals, group differences, and a manifest.

### SCDC (`run_scdc`)

- Engine: R via subprocess, calling `SCDC::SCDC_prop`.
- Required args: `mixture`, `scrna_counts`, `cell_types`, `batch_ids`.
- Optional args: transforms, `ct_sub`, `iter_max`, `nu`, `epsilon`, `weight_basis`, `transform_bisque`, `r_env`, `extra_output_dir`.
- Special features: donor-aware weighted NNLS; SCDC can also do multi-reference ensemble workflows natively, but this wrapper currently exposes the single-reference `SCDC_prop` path.
- Native backend input expectation: bulk and scRNA `Biobase::ExpressionSet` objects with cell type and subject/sample metadata.
- Wrapper preprocessing: validates metadata, normalizes genes, transforms, intersects mixture/scRNA genes, and writes aligned vectors.
- Output: `prop.est.mvw` proportion matrix. Optional extras include `basis.mvw`, `yhat`, `yeval`, and `peval`.

### InstaPrism (`run_instaprism`)

- Engine: R via subprocess, calling `InstaPrism::refPrepare` and `InstaPrism::InstaPrism`.
- Required args: `mixture`, `scrna_counts`, `cell_types`.
- Optional args: transforms, `cell_states`, `r_env`, `extra_output_dir`, `write_z`.
- Special features: fast BayesPrism-like deconvolution with optional cell-state labels; `tool="instaprime"` is accepted as an alias for the original spelling in project notes.
- Native backend input expectation: non-log bulk expression plus scRNA expression with cell type and cell state labels.
- Wrapper preprocessing: validates metadata, normalizes genes, transforms, intersects mixture/scRNA genes, and builds the InstaPrism reference inside R.
- Output: initial posterior cell-type fractions (`t(Post.ini.ct@theta)`). Optional extras can write `theta.tsv` and, if `write_z=True`, `Z_array.rds`.

### TAPE (`run_tape`)

- Engine: Python sidecar subprocess, importing PyTorch-backed `scTAPE`.
- Required args: `mixture`, `scrna_counts`, and `cell_types`.
- Optional args: transforms, `python_env`, `variance_threshold`, scaler, native datatype, gene-length table path, mode/adaptive flags, sparse simulation flag, batch size, epochs, and seed.
- Special features: trains a tissue-adaptive autoencoder from synthetic mixtures generated from the scRNA reference; can also emit adaptive signature output in native overall mode.
- Native backend input expectation: scRNA cells x genes table whose row index is cell type labels, and bulk samples x genes table.
- Wrapper preprocessing: validates metadata, normalizes genes, applies transforms, intersects mixture/scRNA genes, writes TAPE-native orientation, and sets Matplotlib to a noninteractive backend in the sidecar runner.
- Output: sample x cell type fractions from native TAPE prediction. Optional signature output can be written when native TAPE returns a tabular signature matrix.

### BLUE (`run_blue`)

- Engine: Python sidecar subprocess, staging inputs and running upstream BLUE scripts `00_combine_sc_h5ads.py` through `05_predict_bulk.py`.
- Required args: `mixture`, `scrna_counts`, and `cell_types`.
- Optional args: `batch_ids`, transforms, `python_env`, `blue_repo_path`, coarse-to-fine `celltype_mapping`, explicit `input_genes`/`output_genes`, pseudobulk sizes, Dirichlet settings, DEG settings, training hyperparameters, and output paths.
- Special features: predicts cell-type proportions and native BLUE cell-type-specific expression profiles (`predicted_ctGEP_epN.h5ad`). The wrapper forces CPU execution by hiding CUDA devices in the sidecar runner.
- Native backend input expectation: per-library h5ad files with raw counts, a genes x samples bulk TSV whose first column is `gene_id`, a cell-type mapping YAML, and a pipeline YAML.
- Wrapper preprocessing: validates metadata, normalizes gene names, applies transforms, intersects mixture/scRNA genes, writes one h5ad per `batch_ids` group, writes BLUE config files, and uses all shared panel genes as the default BLUE input/output gene lists for Foli-seq compatibility.
- Output: sample x cell type proportion matrix parsed from `predicted_proportions_epN.csv`. Persistent `ctGEP` h5ad output is tracked in `output_paths["ctgep"]` when `ctgep_output_path` or `extra_output_dir` is supplied.

### Scaden (`run_scaden`)

- Engine: Python CLI via subprocess, calling native `scaden simulate`, `scaden process`, `scaden train`, and `scaden predict`.
- Required args: `mixture`; when `retrain=True`, also `scrna_counts` and `cell_types`; when `retrain=False`, `model_dir` is required.
- Optional args: transforms, `python_env`, `model_dir`, `n_training_samples`, `cells_per_sample`, `train_steps`, `batch_size`, `learning_rate`, `var_cutoff`, `seed`, `unknown_celltypes`.
- Special features: trains an ensemble of three neural nets from synthetic pseudobulk generated from the scRNA reference; existing trained models can be reused.
- Native backend input expectation: scRNA cells x genes text plus `Celltype` labels for simulation, mixture TSV genes x samples for prediction, and a trained Scaden model directory.
- Wrapper preprocessing: validates, normalizes gene names, intersects mixture/scRNA genes for training mode, writes Scaden-native temporary files, uses numeric temporary cell indices because native pandas parsing casts the index when `dtype=np.float32`, passes an absolute simulation prefix because native Scaden writes merged h5ad paths relative to the current working directory, and keeps `var_cutoff` default low for panel data.
- Output: sample x cell type proportion matrix from `scaden_predictions.txt`.

### DISSECT (`run_dissect`)

- Engine: Python sidecar subprocess, importing DISSECT and running its simulator, dataset preparation, and fraction-estimation model.
- Required args: `mixture`, `scrna_counts`, `cell_types`.
- Optional args: `batch_ids`, transforms, `python_env`, simulation sizes, training steps, model count, `prop_sparse`, `alpha_range`, native normalization options, `extra_output_dir`.
- Special features: semi-supervised consistency regularization using both synthetic mixtures and real target mixtures; optional batch-aware simulation when `batch_ids` is supplied.
- Native backend input expectation: bulk TSV genes x samples and scRNA h5ad cells x genes with `.obs["cell_type"]` and optional `.obs["batch"]`.
- Wrapper preprocessing: validates metadata, normalizes gene names, applies mixture transform (default `cpm`), intersects mixture/scRNA genes, stages TSV matrix/metadata payloads, then writes h5ad inside the sidecar env so old DISSECT/AnnData can read its own string encoding. The wrapper relaxes default filtering for panel data.
- Output: ensemble `dissect_fractions.txt` as sample x cell type proportions. Optional extras copy `dissect_scores.txt` and DISSECT `main_config.json`.

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
- `DeconvolutionResult.proportion`: `pd.DataFrame | None` with fraction-like sample x output estimates
- `DeconvolutionResult.score`: `pd.DataFrame | None` with non-compositional score outputs
- `DeconvolutionResult.component`: `pd.DataFrame | None` with latent component estimates
- `DeconvolutionResult.p_value`: `pd.DataFrame | None` with p-value outputs when exposed
- `DeconvolutionResult.uncertainty`: `pd.DataFrame | None` with uncertainty outputs when exposed
- `DeconvolutionResult.selected_features`: `pd.DataFrame | None` with selected feature outputs when exposed
- `DeconvolutionResult.output_paths`: `dict[str, str]` containing persistent side-output paths
- `DeconvolutionResult.metadata`: `dict[str, str | int | float | bool]` containing wrapper/runtime settings

Feature-selection wrappers return:

- `FeatureSelectionResult.tool`: `str`
- `FeatureSelectionResult.selected_genes`: `pd.Index`
- `FeatureSelectionResult.ranking`: `pd.DataFrame | None`
- `FeatureSelectionResult.output_paths`: `dict[str, str]` containing persistent side-output paths
- `FeatureSelectionResult.metadata`: `dict[str, str | int | float | bool]` containing wrapper/runtime settings

Output semantics by tool:

- `xCell`: enrichment scores per sample and cell type. These are not compositional fractions and are not constrained to sum to 1.
- `quanTIseq`: fraction-like immune composition output plus `other`; normal outputs sum to 1 per sample.
- `MCP-counter`: abundance scores per population. These are not compositional fractions and are not constrained to sum to 1.
- `PSEA`: marker reference signals per population. These are relative scores scaled by `target_mean`, not compositional fractions.
- `CellCODE`: surrogate proportion variables per sample and cell type. These are latent scores, not compositional fractions, and require bulk sample group labels.
- `ABIS`: robust-regression immune-cell estimates. Native output is percent scale and can be negative because the regression is unconstrained; wrapper output is fraction scale by default.
- `ESTIMATE`: stromal, immune, ESTIMATE, and optional tumor-purity scores. These are not cell-type proportions.
- `AutoGeneS`: `select_features_autogenes` returns selected genes; `run_autogenes` returns feature-selected regression coefficients and selected-feature output. By default the wrapper clips negatives and normalizes rows into proportions.
- `BLADE`: estimated sample-specific fractions from the variational posterior (`ExpF(Beta)`). Optional extras can include best hyperparameters and run IDs.
- `BLUE`: predicted sample-specific proportions are returned in `proportion`; native predicted cell-type-specific expression profiles are written as `ctGEP` h5ad when a persistent output path is provided.
- `CDSeq`: estimated sample-specific proportions. If no reference is supplied, columns are latent components rather than known cell type names; optional extras can include `estGEP`, `estT`, log posterior, and assignment diagnostics.
- `xCell2`: enrichment-like scores per sample and cell type. This is score output, not ratio output; behavior depends on `raw_scores` and spillover options.
- `EPIC`: fraction-like output matrix from native EPIC. Columns can include non-immune/other categories depending on reference mode.
- `dtangle`: estimated mixture fractions from `result$estimates`.
- `DeconRNASeq` wrapper: constrained coefficients from `limSolve::lsei` (non-negative, sum-to-one constraints are enforced in this wrapper).
- `DWLS`: estimated composition matrix from selected solver (`OLS`, `SVR`, or `DampenedWLS`).
- `MuSiC`: weighted fraction estimate (`Est.prop.weighted`) is returned. When `extra_output_dir` is supplied for the SCE/RDS path, supplementary native outputs (`Est.prop.allgene`, `Weight.gene`, `Var.prop`, `r.squared.full`, and a manifest) are also written to disk.
- `MuSiC2`: estimated proportions (`Est.prop`) are returned. Optional extras can include convergence status, iteration count, removed DE genes, and non-converged sample IDs.
- `MEAD`: estimated proportions (`p_hat`) are returned in `proportion`; standard errors (`p_hat_se`) are returned in `uncertainty` when `calc_var=True`. Optional extras can include confidence intervals and group-difference diagnostics.
- `SCDC`: estimated proportions (`prop.est.mvw`) are returned. Optional extras can include basis matrix and expression-fit diagnostics.
- `InstaPrism`: initial posterior cell-type fractions are returned. Optional extras can include deconvolved gene-expression array `Z`, which can be large.
- `TAPE`: neural-network cell type fractions from the native TAPE prediction table; optional signature output is only written for native modes that return a tabular signature.
- `Scaden`: neural-network ensemble cell type fractions from `scaden_predictions.txt`.
- `DISSECT`: neural-network ensemble cell type fractions from `dissect_fractions.txt`; optional extras can include pre-softmax-like score output from `dissect_scores.txt`.
- `MarkerMap`: selected marker genes plus rank/index table. It is a feature selector, not a deconvolution method.
- `scGeneFit`: selected marker genes plus rank/index table. It is a feature selector, not a deconvolution method.
- `Bisque`: only `bulk.props` (deconvolved proportions) is returned.
- `BayesPrism`: returns only `theta` fractions (`get.fraction` output when `update_gibbs=True`, otherwise initial posterior theta). Output can be type-level or state-level depending on `state_or_type`.
- `CIBERSORTx`: returns only cell-type proportion columns. Supplemental columns from native output (`Correlation`, `RMSE`, `P-value`) are removed.

Current limitation:

- Wrappers intentionally normalize primary outputs into explicit fields (`proportion`, `score`, `component`, and optional diagnostics). Method-specific persistent side outputs are tracked in `output_paths` when exposed.

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
| quanTIseq | no explicit wrapper reference intersection because signatures are native | native quanTIseq uses its TIL10 marker/signature resources and removes missing/unusable genes internally | moderate-to-high panel risk if the panel misses many TIL10 genes; smoke run succeeded with full TIL10 signature genes |
| MCP-counter | no explicit wrapper reference intersection because marker tables are native/local | `MCPcounter.estimate` scores available marker genes; native defaults try to download marker tables unless local paths are supplied | score output can run on marker subsets, but panel must include marker genes; provide `genes_path`/`probesets_path` for offline reproducibility |
| PSEA | wrapper normalizes marker gene names but does not build a scRNA signature | `PSEA::marker` uses marker IDs present in the mixture row names to generate population reference signals | suitable only when the Foli panel contains the requested marker genes; output is score-like, not a fraction estimate |
| CellCODE | wrapper intersects `mixture/signature` before native marker tagging | `tagData` further intersects marker source with the bulk reference and `getAllSPVs` uses marker genes with non-zero tag values | score-only; needs enough marker genes per cell type and at least two bulk sample groups |
| ABIS | wrapper normalizes mixture/signature gene names; R layer intersects mixture and ABIS signature genes | Shiny-source logic uses only intersected genes in `MASS::rlm` and has no explicit high-overlap threshold | runnable on panels that contain enough ABIS signature genes; unconstrained regression can produce negative values when overlap/biology is weak |
| ESTIMATE | no explicit custom reference; wrapper writes transformed mixture as GCT | native `estimateScore` computes ssGSEA-like scores over available stromal/immune signature genes and reports overlap in stdout | score-only method; panel assays need enough ESTIMATE signature genes for stable scores |
| AutoGeneS | wrapper intersects `mixture/signature` before native optimization | native `deconvolve` can use available selected genes from the bulk matrix; optimization quality depends on selected-gene search space | panel-compatible when the optimizer is run on the same panel; whole-transcriptome selected genes will not transfer cleanly to a small panel |
| BLADE | wrapper intersects `mixture/signature_mean/signature_sd`, or `mixture/scrna` before building log mean/SD | native BLADE optimizes on provided common genes and can be slow or unstable when marker count is too low relative to cell types | feasible for panels when enough genes remain to estimate both fractions and reference variability; use small hyperparameter grids for smoke tests |
| BLUE | wrapper intersects `mixture/scrna` before staging BLUE h5ad/TSV inputs and writes the shared panel as BLUE's default input/output gene list | native BLUE computes `scRNA ∩ bulk`, generates pseudobulks, and trains on the configured input genes; low overlap weakens signal | panel-compatible when trained on the same panel; do not reuse a whole-transcriptome BLUE checkpoint for panel-only mixtures unless gene lists match |
| CDSeq | optional reference path intersects `mixture/reference` before call; reference-free mode uses mixture genes only | CDSeq expects optional `reference_gep` to have the same number/order of genes as `bulk_data` | reference-free mode avoids mismatch but returns latent components; reference-labeled mode depends on enough shared genes for meaningful assignment |
| xCell2 | wrapper passes `min_shared_genes` (default `0.2`) to `xCell2Analysis` | `xCell2Analysis` intersects with `getGenesUsed(xcell2object)`; stops if shared fraction `< minSharedGenes`, warns if `< 0.85` | high risk if reference is whole-transcriptome and mixture is panel; prefer panel-matched xCell2 object |
| EPIC | custom-reference path intersects `mixture/signature` before call | `EPIC::EPIC` intersects bulk/reference genes internally; warns when common genes are low (`<2000`); stops if signature genes in common are fewer than reference cell types | moderate risk with very small overlap or too many cell types |
| dtangle | wrapper intersects `mixture/signature` before call | marker selection runs on provided genes; no hard overlap threshold in `dtangle::dtangle` itself | moderate risk; low overlap means weaker marker set |
| DeconRNASeq wrapper (`limSolve::lsei`) | wrapper intersects `mixture/signature` before call | constrained least-squares solve; no explicit overlap threshold | moderate risk when genes are too few relative to number of cell types |
| DWLS | wrapper intersects `mixture/signature`; R layer intersects again | solver runs on intersected genes; no explicit high-overlap threshold in wrapper | moderate risk; marker-poor panels can destabilize estimates |
| MuSiC | matrix path intersects `mixture/scrna` before call; SCE/RDS path keeps the full reference and lets MuSiC intersect internally | `MuSiC::music_prop` checks common genes and can stop with `Too few common genes!`; full-reference SCE/RDS mode can produce different estimates than pre-subset panel references | usually tolerant, but for exact replication prefer full-reference SCE/RDS mode with the MuSiC R 4.4 sidecar |
| MuSiC2 | matrix path intersects `control/case/scrna`; SCE/RDS path intersects only control/case before call | `MuSiC2` uses genes shared by bulk and scRNA, and can remove DE genes iteratively | panel risk depends on enough shared genes remaining after DE-gene filtering |
| MEAD | wrapper intersects `mixture/scrna` before building the reference SCE | native MEAD can filter genes by expression threshold/quantiles; wrapper defaults `filter_gene=False` and `gene_thresh=0.0` for panel compatibility | panel-compatible when enough informative genes remain; re-enable native filtering only when the panel is large enough |
| SCDC | wrapper intersects `mixture/scrna` before call | `SCDC_prop` stops with `Too few common genes!` when common genes are below 20% of the smaller input gene dimension | can be strict for tiny panels; best used with panel-aware references or enough marker coverage |
| InstaPrism | wrapper intersects `mixture/scrna` before building `refPrepare` reference | reference preparation and deconvolution run only on provided shared genes; no wrapper-side high-overlap threshold | likely runnable on panels, but BayesPrism-like performance depends strongly on marker coverage and cell-state labels |
| TAPE | wrapper intersects `mixture/scrna` before native simulation/training | native TAPE performs additional variance filtering independently on simulated and target data before intersecting genes | panel-compatible only if enough informative genes survive `variance_threshold`; for small panels keep `variance_threshold` high and use small `epochs` for smoke tests |
| Scaden | training mode intersects `mixture/scrna` before native simulation; prediction-only mode relies on the trained model gene list | native `process` intersects training/prediction genes and filters by variance; native `predict` requires all model genes in the mixture | panel-compatible if the model is trained/processed against the same panel; prediction-only with a whole-transcriptome model is likely to fail |
| DISSECT | wrapper intersects `mixture/scrna` before h5ad/TSV creation | native `prepare_data` intersects simulated and target genes again; native filtering can remove genes if `var_cutoff`/QC filters are strict | panel-compatible in principle, but training is heavy and performance depends on enough informative genes after native preprocessing |
| Bisque | wrapper intersects `mixture/scrna` before call | `BisqueRNA` overlap utilities require non-zero overlap and can stop if filtering leaves zero genes | moderate risk if filtering removes most panel genes |
| BayesPrism | wrapper intersects `mixture/scrna` before call | `BayesPrism::new.prism` intersects again; stops at zero shared genes; warns when shared genes `< 100` | moderate risk for very tiny panels; ~500-gene panels generally satisfy this warning threshold |
| CIBERSORTx | wrapper intersects `mixture/signature` before container run | native fractions image overlap minimum is not explicitly surfaced by wrapper; output computed on intersected genes | unknown-to-moderate risk; empirical wrapper runs succeeded with 20 and 120 shared genes in this repo |

Why this matters for Foli-seq:
- panel assays often keep only a small subset of transcriptome genes
- methods with hard overlap requirements (`xCell`, strict `xCell2` settings) are most likely to fail
- for other methods, runs may complete but precision can degrade as overlap/marker informativeness drops

## Backends

- R-backed wrappers call `src/foli_decon/r_scripts/run_tools.R` via `Rscript` subprocess
- TAPE, BLADE, BLUE, AutoGeneS, Scaden, and DISSECT wrappers call Python sidecar commands via subprocess
- CIBERSORTx wrapper calls containerized `cibersortx/fractions` via `podman`/`docker` subprocess
- current implementation does not use `rpy2`; R communication is file + subprocess based
- shared preprocess logic lives in `src/foli_decon/preprocess.py`

## Notes on tool-specific extras

- `xCell2` can require lineage-specific references (for example immuno-oncology); pass the correct `.rds` object.
- `MCP-counter` may try to download marker tables through its native defaults; pass `genes_path` or `probesets_path` for offline/reproducible runs.
- `PSEA` requires explicit marker sets and returns relative reference signals rather than fractions.
- `CellCODE` requires bulk sample group labels and returns SPV scores rather than fractions.
- `AutoGeneS` is best run in its sidecar env because its conda dependency set includes old AnnData-era packages.
- `BLADE` is best run in its sidecar env because the native package is `BLADE-Deconvolution`, not the unrelated PyPI package named `BLADE`.
- `BLUE` is best run in its CPU sidecar because upstream defaults use `uv` with CUDA-oriented Torch wheels, while this repo needs CPU-only execution by default.
- `CDSeq` can run without a reference, but then output columns are latent components, not known cell type labels.
- `MuSiC` and `Bisque` require donor/sample IDs (`batch_ids`) for reference cells.
- `MuSiC2` requires separate control and case bulk matrices; it is not part of the default single-mixture two-split benchmark CLI.
- `SCDC` requires donor/sample IDs and can be strict about common-gene fraction.
- `InstaPrism` is best treated as a sidecar method because it is GitHub-only and compiles Rcpp code.
- `TAPE` is best run in its CPU sidecar so conda supplies `pytorch-cpu` and pip does not pull CUDA Torch wheels.
- `Scaden` and `DISSECT` are CPU TensorFlow sidecar methods; start with small `train_steps`/`n_training_samples` for smoke tests.
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

print(result.proportion.head())
```

Notes:
- CIBERSORTx container runs require `username` and `token` (password is not used by the fractions image).
- The wrapper aligns mixture/signature genes before writing container input files.

## Current strategy for dependency conflicts

Primary path:
1. single environment `foli-decon`
2. conda/mamba packages first for all available decon tools
3. small post-install script only for packages missing from conda (currently `xCell2`)
4. affiliated sidecar envs for incompatible R, TensorFlow, PyTorch, or old Python stacks (`SCDC`, `InstaPrism`, `CDSeq`, `TAPE`, `BLADE`, `BLUE`, `Scaden`, `DISSECT`, `AutoGeneS`)

Fallback path if a method remains incompatible in the unified env:
1. keep Python wrapper API stable
2. run that method in a dedicated container
3. return the same `samples x cell_types` output format
