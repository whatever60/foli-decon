# Project Contract And TODO

This file preserves decisions from the working thread so future work does not
lose direction after context compaction. Treat this as the project contract
unless the user explicitly changes it.

## Core Goal

`foli-decon` wraps cell type deconvolution tools for Foli-seq panel RNA-seq.
Foli-seq inputs are UMI count matrices over a targeted panel of roughly 500
genes. Gene mismatch with scRNA-seq or bulk transcriptome references is normal,
not exceptional.

The package should provide a Python-facing API that hides backend details while
remaining explicit about each method's native expectations, transformations,
extra metadata, gene-overlap behavior, and output semantics.

## API Contract

- Use singular result fields: `result.proportion`, `result.score`,
  `result.component`, `result.p_value`, `result.uncertainty`, and
  `result.selected_features`.
- Do not keep a `.proportions` compatibility alias. This branch intentionally
  made a hard break.
- Missing outputs are `None`, not empty data frames.
- Methods that return scores should populate `score`, not `proportion`.
- Built-in score/signature methods should not expose `train_reference`.
- Reference-generation APIs should exist only when the method has a real native
  reference-building concept.
- Feature selectors should return `FeatureSelectionResult`; end-to-end helpers
  may return two objects, but should not create one busy wrapper dataclass.
- Prefer explicit arguments over opaque folders such as `reference_dir`.

## Backend Contract

- Do not directly import or call `immunedeconv` or `omnideconv`.
- Their source can be read for reference, but wrappers should use native
  package APIs whenever practical.
- R tools may use R subprocesses and small R helper scripts.
- CIBERSORTx should stay containerized through `podman`/`docker`.
- Python tools with incompatible stacks should use sidecar subprocesses.

## Environment Contract

- Work on `add-more-decon-tools`; do not touch `main`.
- Do not mutate existing conda/mamba envs during expansion/probe work unless
  the user explicitly allows it.
- Prefer conda/mamba packages over source installs.
- Source installs are acceptable inside sidecar envs when dependencies are
  mostly pinned by conda first.
- Combine tools into as few envs as is sane, but prefer stable sidecars over
  forcing ancient or incompatible dependencies into the main env.
- CPU-only is the default for deep-learning tools.
- Use temporary/probe envs only when dry-runs are insufficient to inspect native
  APIs.

## Gene-Mismatch Contract

- Prefer method-native gene matching/filtering when the tool clearly provides
  it.
- Use gene intersection as the fallback.
- Do not add generic `min_shared_genes` or `min_shared_fraction` controls unless
  the method already exposes that concept or the wrapper already has that knob.
- If behavior is ambiguous, fail loudly.

## Testing Contract

- Keep normal pytest tests small and synthetic.
- Optional smoke scripts may use local data under `~/dev`.
- Small real-data fixtures can be copied into the repo only when they materially
  improve regression coverage.

## Method Decisions

- `TIMER`: conceptually kept, but local native support is unresolved and is
  explicitly deferred.
- `TIMER2` and `TIMER3`: web-based routes are not local wrappers and should stay
  unsupported.
- `AutoGeneS`: primarily a feature/gene selection method feeding downstream
  deconvolution; `run_autogenes` may expose its native deconvolution workflow,
  but feature selection must be available as a first-class API.
- `PSEA`: score output, not compositional proportions.
- `TAPE`: PyPI package is `scTAPE`; use the CPU PyTorch sidecar.
- `Scaden`: supported through the shared CPU TensorFlow sidecar. Wrapper must
  use numeric temporary cell indices and an absolute simulation prefix to avoid
  Scaden 1.1.2 native I/O quirks.
- `DISSECT`: supported through the shared CPU TensorFlow sidecar. Stage scRNA
  matrix/metadata as TSV and let the sidecar runner write h5ad with its own
  old AnnData version.
- `BLADE`: PyPI package is `BLADE-Deconvolution`; do not use the unrelated
  package named `BLADE`.
- `MEAD`: supported through a sidecar R env plus GitHub source installs for
  `vashr` and `MEAD`. It is scRNA-reference based, returns `proportion`, and
  can also return `uncertainty`.
- `CellCODE`: supported through a sidecar R env plus GitHub source install.
  It returns surrogate proportion variable scores in `score`, not fractions.
- `ABIS`: supported through local source logic copied from the ABIS Shiny app
  pattern, using local ABIS signature resources and `MASS::rlm`. Do not fetch
  the signatures at runtime.
- `ESTIMATE`: supported through the original R-Forge `estimate` package. It is
  a source install, but it is no-compile and narrow enough for the post-install
  script.
- `BLUE`: supported through a CPU sidecar pipeline adapter. The wrapper stages
  h5ad/TSV/YAML inputs, runs upstream BLUE scripts `00` through `05`, returns
  proportions, and records persistent ctGEP h5ad output paths when requested.
- `MarkerMap`: feature selector only. Supported through
  `environment-py311-cpu.yml` plus pip/source install layer.
- `scGeneFit`: feature selector only. Supported through the same feature
  selector sidecar.
- `DSA`: unresolved/deferred until a stable local package/source route is
  identified.

## Current TODO

1. Keep the current branch contract-aligned after each new wrapper.
2. Convert the ad hoc sidecar smoke commands for MarkerMap, scGeneFit, TAPE,
   BLUE, Scaden, and DISSECT into checked-in optional smoke scripts if we want
   repeatable local validation.
3. Continue source/package triage for unresolved methods one by one.
4. Keep benchmark defaults conservative; sidecar-heavy, credentialed,
   marker-set, or non-single-mixture tools should be explicit opt-ins.
5. Preserve NEC bit-identity constraints in `/home/ubuntu/dev/20251201_nec`
   work; only package-version differences are acceptable explanations for
   non-identical results.
