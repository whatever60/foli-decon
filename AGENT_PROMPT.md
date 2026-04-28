# AGENT_PROMPT

This document summarizes the project logic and execution roadmap as an agent prompt that can be reused by collaborators.

## Copy/Paste Agent Prompt

You are implementing `foli-decon`, a Python package for cell-type deconvolution on Foli-seq panel RNA-seq data.

### Project purpose

Build a unified interface to benchmark multiple existing deconvolution tools on Foli-seq UMI count matrices using scRNA-seq or cell-type-specific bulk references.

Foli-seq assumptions:
- panel RNA-seq (~500 genes)
- UMI count input
- frequent gene mismatch between mixture and reference

### Target tools

Implement and maintain wrappers for:
- xCell
- xCell2
- EPIC
- dtangle
- BayesPrism
- DeconRNASeq
- CIBERSORTx
- Bisque
- DWLS
- MuSiC

### Non-negotiable interface goals

1. Keep one consistent Python API for all tools (`tool="..."` dispatcher + per-tool wrappers).
2. Accept `pandas` matrices/vectors in memory; wrapper handles file conversion for backend engines.
3. Automatically handle gene mismatch by intersection before backend execution.
4. Clearly document required and optional tool-specific inputs (for example `batch_ids`, `cell_states`, xCell2 object path, CIBERSORTx credentials).
5. Keep wrapper-side transformations explicit and consistent (`raw`, `cpm`, `cpm_log1p`, `tpm`, `tpm_log1p`).
6. Standardize primary output into `samples x outputs` DataFrame inside a `DeconvolutionResult` object.
7. Explicitly document output semantics: score vs proportion, and what native supplemental statistics are dropped or retained.

### Backend strategy

- Prefer conda/mamba installation for heavy R dependencies first.
- Minimize source compilation from R when possible.
- Use subprocess-based R execution (`Rscript`) with temporary TSV/JSON payload files.
- Use container runtime (`podman`/`docker`) for CIBERSORTx.
- If unified env becomes impractical, keep Python API stable and isolate problematic tools in per-tool containers.

### External code references to reuse

Use existing local examples before writing new complex code:
- `~/dev/20251201_nec` (xCell2, MuSiC, plotting workflows)
- `~/dev/exfoseq_deconvolution` (CIBERSORTx via podman and input construction)

### Engineering principles

- Keep design simple and modular.
- Reuse existing code patterns across wrappers.
- Fail loudly on invalid inputs rather than silently correcting.
- Keep preprocessing and I/O behavior centralized in shared utilities.
- Avoid backend-specific surprises in high-level API.

### Phased roadmap

Phase 1: Environment and dependencies
- Maintain conda/mamba env `foli-decon`.
- Pin versions to keep R stack compatible.
- Install non-conda stragglers in a minimal post-install step.

Phase 2: Shared core
- Define shared preprocessing (validation, gene normalization, transforms, intersections).
- Define shared result model (`DeconvolutionResult`).
- Build unified dispatcher and per-tool entrypoints.

Phase 3: Wrapper implementation
- Implement R-backed wrappers (xCell, xCell2, EPIC, dtangle, DeconRNASeq, DWLS, MuSiC, Bisque, BayesPrism).
- Implement CIBERSORTx container wrapper.
- Normalize outputs to the common return object.

Phase 4: Documentation
- For each tool document:
  - execution engine (Python direct, R subprocess, container subprocess)
  - arguments (required/optional, expected content)
  - special capabilities (for example donor IDs, cell states)
  - native input expectations and wrapper-applied transforms
  - output semantics and supplemental metrics behavior

Phase 5: Validation
- Unit-test wrapper plumbing and preprocessing behavior.
- Mock external engines in unit tests where needed.
- Perform at least one real CIBERSORTx run to validate true container path.

Phase 6: Fallback hardening
- For any tool that cannot coexist in unified env, isolate in container and preserve API contract.

### Acceptance criteria

- All listed tools have callable wrappers under one Python package.
- Gene intersection behavior is consistent and documented.
- Tool-specific required metadata is enforced (`batch_ids`, `cell_types`, etc.).
- Output contract is consistent, with clear documentation for exceptions (xCell/xCell2 scores).
- Environment setup steps are reproducible.
- Tests pass in `foli-decon` environment.

## Snapshot of what is already implemented

Current repo status already includes:
- unified dispatcher and per-tool wrappers for all 10 target tools
- shared preprocessing/transforms/intersection utilities
- subprocess-based R backend runner
- CIBERSORTx container wrapper and real end-to-end run validation
- README API documentation with per-tool engine/input/output notes
- passing test suite in `foli-decon`
