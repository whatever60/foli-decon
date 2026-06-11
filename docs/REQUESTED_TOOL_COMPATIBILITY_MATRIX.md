# Requested Tool Compatibility Matrix

Last updated: 2026-06-11 on branch `add-more-decon-tools`.

This is the explicit compatibility matrix promised during the method-expansion
plan. It consolidates the package probes, dry-run decisions, sidecar decisions,
wrapper status, and known blockers for the methods requested beyond the original
core wrapper set.

Status meanings:

- `supported`: Python-facing wrapper exists in this branch.
- `feature selector`: wrapper returns selected genes, not deconvolution.
- `deferred`: route was investigated, but no local wrapper is exposed yet.
- `excluded`: intentionally not pursued in this branch.

## Supported Or Integrated Methods

| Method | Role and output | Package/source route checked | Env decision | Wrapper/API status | Validation status | Main compatibility notes |
|---|---|---|---|---|---|---|
| xCell | Built-in signature score | conda `r-xcell` | main `foli-decon` R 4.2 env | `run_xcell`, returns `score` | unit tests | Native signatures only; panel data must contain enough marker genes for meaningful scores. |
| xCell2 | Trained signature score | GitHub `dviraran/xCell2`; conda dependencies for heavy imports | main env plus `scripts/install_r_packages.R` | `train_xcell2_reference`, `run_xcell2`, returns `score` | unit tests | Requires trained xCell2 object or scRNA training input; lineage/cell ontology inputs stay explicit when used. |
| MCP-counter | Built-in marker score | conda `r-mcpcounter` | main env | `run_mcp_counter`, returns `score` | unit tests with local marker tables; `r-base=4.2 r-mcpcounter` dry-run solved | Native defaults can download marker tables at runtime; use local `genes_path` or `probesets_path` for reproducibility. |
| quanTIseq | Built-in immune fractions | conda `bioconductor-quantiseqr` | main env | `run_quantiseq`, returns `proportion` | unit tests and small R smoke noted in README history | Built-in immune reference; no custom scRNA reference training. |
| ConsensusTME | Built-in TME score | conda `r-consensustme` | main env | `run_consensus_tme`, returns `score` | unit tests | Score method, not compositional proportions. |
| PSEA | Marker reference-signal score | conda `bioconductor-psea` | main env | `run_psea`, returns `score` | unit tests | Requires user marker sets; not a fraction estimator. |
| ABIS | Signature-based fractions | ABIS Shiny source/resources; `MASS::rlm` logic | main env plus local ABIS resources | `run_abis`, returns `proportion` | unit tests | No conda package; run `scripts/download_abis_resources.py` or pass `signature_path`. Native coefficients are percent and are rescaled to fractions by default. |
| ESTIMATE | Stromal/immune/tumor purity scores | R-Forge `estimate 1.0.13`; no compile | main env post-install through `scripts/install_r_packages.R` | `run_estimate`, returns `score` | unit tests | Source install is acceptable because it is narrow and no-compile. |
| CellCODE | Surrogate proportion variable scores | GitHub `mchikina/CellCODE`; conda `gplots` and `sva` | `environment-r44.yml` | `run_cellcode`, returns `score` | unit tests; merged R44 env created; package load passed | Native output is SPV-like scores, not cell fractions. |
| AutoGeneS | Feature selection plus optional native deconvolution | conda `autogenes=1.0.4` | `environment-py311-cpu.yml` | `select_features_autogenes` returns selected genes; `run_autogenes` returns deconvolution output | unit tests; merged Python env created; import passed | Kept out of main env because conda solve pulls old AnnData-era dependencies. |
| MarkerMap | Feature selector | PyPI `markermap`; GitHub `iancovert/persist`; conda PyTorch/Lightning/Scanpy stack | `environment-py311-cpu.yml` | `select_features_markermap`, returns `FeatureSelectionResult` | real synthetic smoke test selected markers in merged env | Feature selector only; no deconvolution layer is exposed. |
| scGeneFit | Feature selector | PyPI `scGeneFit==1.0.2` | `environment-py311-cpu.yml` | `select_features_scgenefit`, returns `FeatureSelectionResult` | real tiny synthetic smoke test selected markers in merged env | Installed with conda-provided NumPy/SciPy/scikit-learn to avoid deprecated PyPI dependency behavior. |
| BLADE | Signature/scRNA deconvolution fractions | PyPI `BLADE-Deconvolution==0.0.7`; PyPI `BLADE` is unrelated | `environment-py311-cpu.yml` | `run_blade`, returns `proportion` | unit tests; real synthetic smoke test passed in merged env | Python sidecar avoids forcing package-specific numerical pins into main env. |
| BLUE | Deep-learning deconvolution fractions plus ctGEP h5ad | GitHub `SichenZhu/BLUE`; conda CPU PyTorch/Scanpy stack | `environment-py311-cpu.yml` | `run_blue`, returns `proportion`; persistent ctGEP h5ad is tracked in `output_paths["ctgep"]` | unit tests; real tiny CPU smoke test completed with `epochs=1` in merged env | Wrapper stages Foli-decon inputs into BLUE h5ad/TSV/YAML pipeline inputs, runs native scripts `00` through `05`, and forces CPU execution by hiding CUDA devices. |
| TAPE | Deep-learning deconvolution fractions | PyPI `scTAPE==1.1.2`; no R/conda alias | `environment-py311-cpu.yml` | `run_tape`, returns `proportion` | real tiny CPU smoke test with `epochs=1` passed in merged env | CPU-only PyTorch sidecar; runner restores `torch.load(weights_only=False)` for native TAPE's same-run checkpoint under PyTorch 2.6+. |
| Scaden | Deep-learning deconvolution fractions | PyPI `scaden==1.1.2` | `environment-deep-cpu.yml` | `run_scaden`, returns `proportion` | real tiny CPU smoke test with `train_steps=1` | Shared TensorFlow 2.7 sidecar; wrapper handles native numeric cell-ID and absolute-prefix I/O quirks. |
| DISSECT | Deep-learning deconvolution fractions | GitHub `imsb-uke/DISSECT` | `environment-deep-cpu.yml` | `run_dissect`, returns `proportion`; optional scores copied to extras | real tiny CPU smoke test with `train_steps=1` | Shared TensorFlow sidecar; wrapper stages TSV and lets sidecar write old-AnnData-compatible h5ad. |
| CDSeq | Reference-free or reference-assisted fractions/components | conda `r-cdseq` | `environment-r44.yml` | `run_cdseq`, returns `proportion` or `component` depending native mode | unit tests; merged R44 package load passed | Kept outside main env because available conda build targets R 4.4 and pulls a large Seurat/Harmony stack. |
| MuSiC | scRNA-reference fractions | conda `r-music` | main env; optional R 4.4 sidecar for exact MuSiC 1.0.0 workflows | `run_music`, returns `proportion` | unit tests and prior NEC workflow checks; merged R44 package load passed | Main env stays R 4.2 for BayesPrism; exact newer MuSiC workflows use `environment-r44.yml`. |
| MuSiC2 | Case/control scRNA-reference differential deconvolution | GitHub `xuranw/MuSiC` MuSiC 1.0.0 route | `environment-r44.yml` recommended | `run_music2`, returns `proportion` plus optional extras | unit tests; merged R44 package load passed | Requires separate control and case mixtures; not included in generic single-mixture benchmark defaults. |
| MEAD | scRNA-reference fractions plus uncertainty | GitHub `DongyueXie/MEAD`; GitHub dependency `mengyin/vashr` | `environment-r44.yml` | `run_mead`, returns `proportion` and `uncertainty` when native output provides it | unit tests; merged R44 package load passed | No conda MEAD/vashr packages; compiled dependencies are conda-managed, source installs are limited to MEAD/vashr. |
| SCDC | scRNA-reference ensemble fractions | conda `r-scdc` | `environment-scdc-r40.yml` | `run_scdc`, returns `proportion` | unit tests; sidecar dry-run solved | Isolated because available conda build requires R 4.0 and conflicts with the main R 4.2/BayesPrism stack. |
| InstaPrism | scRNA-reference fractions | GitHub InstaPrism package route | `environment-r44.yml` | `run_instaprism`, returns `proportion` | unit tests; merged R44 package load passed | User alias `instaprime` maps to `instaprism`; source install is isolated in sidecar. |

## Deferred Or Excluded Methods

| Method | Probe result | Decision | Reason |
|---|---|---|---|
| DSA | `r-dsa` absent; CRAN/PyPI names found were unrelated | deferred | No stable local deconvolution package/source route was identified. |
| TIMER | `r-timer` absent | deferred | Local native implementation route unresolved. Web-only routes are not acceptable for this package interface. |
| TIMER2 | web tool | excluded | Web service route only; not a local package API. |
| TIMER3 | web tool / no local package | excluded | Web service route only; removed from supported-tool semantics. |
| Linseed | user later chose to skip | excluded | Not part of this branch's target set. |
| GS-NMF | user later chose to skip | excluded | Not part of this branch's target set. |

## BLUE Adapter Details

BLUE is now treated as a supported pipeline-backed method rather than deferred.
The method trains on scRNA-derived pseudobulks and predicts both cell-type
proportions and cell-type-specific expression profiles from bulk input.

The wrapper intentionally wraps multiple native source commands in one
Foli-decon API call. `run_blue` performs the following adapter work:

- stage one or more per-donor/per-sample scRNA h5ad files with raw counts in
  `.X` or a configured raw-count layer
- stage a genes-by-samples bulk TSV with first column `gene_id`
- generate a BLUE `celltype_mapping.yaml`
- generate a BLUE `pipeline_config.yaml` with all paths rewritten into a temp
  or persistent output directory
- run `scripts/deconv/00_combine_sc_h5ads.py` through
  `scripts/deconv/05_predict_bulk.py` in order
- optionally run `06_visualize_predictions.py` or
  `06_eval_held_out_pseudobulks.py` only for diagnostics, not core output
- collect `predicted_proportions_epN.csv` into `DeconvolutionResult.proportion`
- collect `predicted_ctGEP_epN.h5ad` into `output_paths` because this is a
  native BLUE-specific extra output that most other tools do not produce
- force CPU-safe Torch installation in a sidecar instead of using upstream
  CUDA-12.8 uv defaults
- smoke-test the whole path on a tiny synthetic scRNA/bulk fixture

Current validation status: unit coverage exists for staging, payload generation,
dispatcher support, and output parsing. A real CPU sidecar smoke test completed
on a tiny synthetic scRNA/bulk fixture with `epochs=1`; it produced
`cache/blue_smoke/blue_proportions.tsv` and
`cache/blue_smoke/blue_predicted_ctGEP.h5ad`.

## Reproducibility Commands

The package/source availability sweep can be repeated with:

```bash
scripts/check_requested_method_availability.sh
```

The sidecar env specs created for supported-but-incompatible methods are:

```text
environment-r44.yml
environment-deep-cpu.yml
environment-py311-cpu.yml
environment-scdc-r40.yml
```

Legacy per-tool sidecar YAMLs are retained as troubleshooting/probe artifacts,
but the compressed envs above are the recommended installation targets.
