"""Wrappers for deconvolution tools that directly use scRNA-seq references."""

import tempfile
from pathlib import Path

import pandas as pd

from foli_decon.io import write_expression_matrix, write_vector
from foli_decon.models import DeconvolutionResult
from foli_decon.preprocess import (
    align_mixture_and_scrna,
    apply_transform,
    normalize_gene_names,
    validate_expression_matrix,
)
from foli_decon.r_runner import run_r_tool


def _validate_scrna_metadata(
    scrna_counts: pd.DataFrame,
    cell_types: pd.Series,
    batch_ids: pd.Series,
) -> None:
    """Validate that scRNA metadata vectors align to cell columns."""

    if len(cell_types) != scrna_counts.shape[1]:
        raise ValueError("cell_types length must match number of scRNA columns")
    if len(batch_ids) != scrna_counts.shape[1]:
        raise ValueError("batch_ids length must match number of scRNA columns")


def run_music(
    mixture: pd.DataFrame,
    scrna_counts: pd.DataFrame | None = None,
    cell_types: pd.Series | None = None,
    batch_ids: pd.Series | None = None,
    mixture_transform: str = "raw",
    scrna_transform: str = "raw",
    select_ct: list[str] | None = None,
    scrna_sce_rds_path: str | None = None,
    cell_type_column: str = "cellType",
    batch_id_column: str = "sampleID",
    r_env: str | None = None,
    r_output_path: str | None = None,
    extra_output_dir: str | None = None,
) -> DeconvolutionResult:
    """Run MuSiC using bulk mixtures plus single-cell reference counts."""

    validate_expression_matrix(mixture, "mixture")
    rscript_command = None
    if r_env is not None:
        rscript_command = ["mamba", "run", "-n", r_env, "Rscript"]

    if scrna_sce_rds_path is not None:
        if scrna_transform != "raw":
            raise ValueError("scrna_transform must be 'raw' when scrna_sce_rds_path is used")
        mixture_norm = normalize_gene_names(mixture)
        mixture_transformed = apply_transform(mixture_norm, transform=mixture_transform)

        with tempfile.TemporaryDirectory(prefix="foli_music_sce_") as workdir:
            workdir_path = Path(workdir)
            mixture_path = workdir_path / "mixture.tsv"
            write_expression_matrix(mixture_transformed, str(mixture_path))
            result, _ = run_r_tool(
                tool="music_sce",
                payload={
                    "mixture_path": str(mixture_path),
                    "scrna_sce_rds_path": scrna_sce_rds_path,
                    "cell_type_column": cell_type_column,
                    "batch_id_column": batch_id_column,
                    "select_ct": select_ct,
                    "extra_output_dir": extra_output_dir,
                },
                rscript_command=rscript_command,
                output_path=r_output_path,
            )

        return DeconvolutionResult(
            tool="music",
            proportions=result,
            metadata={
                "mixture_transform": mixture_transform,
                "scrna_transform": scrna_transform,
                "select_ct_count": 0 if select_ct is None else len(select_ct),
                "reference_mode": "sce_rds",
                "cell_type_column": cell_type_column,
                "batch_id_column": batch_id_column,
                "r_env": "" if r_env is None else r_env,
                "r_output_path": "" if r_output_path is None else r_output_path,
                "extra_output_dir": "" if extra_output_dir is None else extra_output_dir,
            },
        )

    if scrna_counts is None or cell_types is None or batch_ids is None:
        raise ValueError("Provide scrna_counts, cell_types, and batch_ids, or provide scrna_sce_rds_path")

    validate_expression_matrix(scrna_counts, "scrna_counts")
    _validate_scrna_metadata(scrna_counts, cell_types, batch_ids)

    mixture_norm = normalize_gene_names(mixture)
    scrna_norm = normalize_gene_names(scrna_counts)
    mixture_aligned, scrna_aligned = align_mixture_and_scrna(mixture_norm, scrna_norm)

    mixture_transformed = apply_transform(mixture_aligned, transform=mixture_transform)
    scrna_transformed = apply_transform(scrna_aligned, transform=scrna_transform)

    with tempfile.TemporaryDirectory(prefix="foli_music_") as workdir:
        workdir_path = Path(workdir)
        mixture_path = workdir_path / "mixture.tsv"
        scrna_path = workdir_path / "scrna.tsv"
        cell_types_path = workdir_path / "cell_types.tsv"
        batch_ids_path = workdir_path / "batch_ids.tsv"

        write_expression_matrix(mixture_transformed, str(mixture_path))
        write_expression_matrix(scrna_transformed, str(scrna_path))
        write_vector(pd.Series(cell_types.to_numpy()), str(cell_types_path), "cell_type")
        write_vector(pd.Series(batch_ids.to_numpy()), str(batch_ids_path), "batch_id")

        result, _ = run_r_tool(
            tool="music",
            payload={
                "mixture_path": str(mixture_path),
                "scrna_path": str(scrna_path),
                "cell_types_path": str(cell_types_path),
                "batch_ids_path": str(batch_ids_path),
                "select_ct": select_ct,
            },
            rscript_command=rscript_command,
            output_path=r_output_path,
        )

    return DeconvolutionResult(
        tool="music",
        proportions=result,
        metadata={
            "mixture_transform": mixture_transform,
            "scrna_transform": scrna_transform,
            "select_ct_count": 0 if select_ct is None else len(select_ct),
            "reference_mode": "matrix",
            "r_env": "" if r_env is None else r_env,
            "r_output_path": "" if r_output_path is None else r_output_path,
        },
    )


def run_bisque(
    mixture: pd.DataFrame,
    scrna_counts: pd.DataFrame,
    cell_types: pd.Series,
    batch_ids: pd.Series,
    mixture_transform: str = "raw",
    scrna_transform: str = "raw",
    use_overlap: bool = False,
    old_cpm: bool = True,
) -> DeconvolutionResult:
    """Run Bisque RNA ReferenceBasedDecomposition."""

    validate_expression_matrix(mixture, "mixture")
    validate_expression_matrix(scrna_counts, "scrna_counts")
    _validate_scrna_metadata(scrna_counts, cell_types, batch_ids)

    mixture_norm = normalize_gene_names(mixture)
    scrna_norm = normalize_gene_names(scrna_counts)
    mixture_aligned, scrna_aligned = align_mixture_and_scrna(mixture_norm, scrna_norm)

    mixture_transformed = apply_transform(mixture_aligned, transform=mixture_transform)
    scrna_transformed = apply_transform(scrna_aligned, transform=scrna_transform)

    with tempfile.TemporaryDirectory(prefix="foli_bisque_") as workdir:
        workdir_path = Path(workdir)
        mixture_path = workdir_path / "mixture.tsv"
        scrna_path = workdir_path / "scrna.tsv"
        cell_types_path = workdir_path / "cell_types.tsv"
        batch_ids_path = workdir_path / "batch_ids.tsv"

        write_expression_matrix(mixture_transformed, str(mixture_path))
        write_expression_matrix(scrna_transformed, str(scrna_path))
        write_vector(pd.Series(cell_types.to_numpy()), str(cell_types_path), "cell_type")
        write_vector(pd.Series(batch_ids.to_numpy()), str(batch_ids_path), "batch_id")

        result, _ = run_r_tool(
            tool="bisque",
            payload={
                "mixture_path": str(mixture_path),
                "scrna_path": str(scrna_path),
                "cell_types_path": str(cell_types_path),
                "batch_ids_path": str(batch_ids_path),
                "use_overlap": use_overlap,
                "old_cpm": old_cpm,
            },
        )

    return DeconvolutionResult(
        tool="bisque",
        proportions=result,
        metadata={
            "mixture_transform": mixture_transform,
            "scrna_transform": scrna_transform,
            "use_overlap": use_overlap,
            "old_cpm": old_cpm,
        },
    )


def run_bayesprism(
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
) -> DeconvolutionResult:
    """Run BayesPrism deconvolution with single-cell references."""

    validate_expression_matrix(mixture, "mixture")
    validate_expression_matrix(scrna_counts, "scrna_counts")
    if len(cell_types) != scrna_counts.shape[1]:
        raise ValueError("cell_types length must match number of scRNA columns")
    if cell_states is not None and len(cell_states) != scrna_counts.shape[1]:
        raise ValueError("cell_states length must match number of scRNA columns")

    mixture_norm = normalize_gene_names(mixture)
    scrna_norm = normalize_gene_names(scrna_counts)
    mixture_aligned, scrna_aligned = align_mixture_and_scrna(mixture_norm, scrna_norm)

    mixture_transformed = apply_transform(mixture_aligned, transform=mixture_transform)
    scrna_transformed = apply_transform(scrna_aligned, transform=scrna_transform)

    with tempfile.TemporaryDirectory(prefix="foli_bayesprism_") as workdir:
        workdir_path = Path(workdir)
        mixture_path = workdir_path / "mixture.tsv"
        scrna_path = workdir_path / "scrna.tsv"
        cell_types_path = workdir_path / "cell_types.tsv"

        write_expression_matrix(mixture_transformed, str(mixture_path))
        write_expression_matrix(scrna_transformed, str(scrna_path))
        write_vector(pd.Series(cell_types.to_numpy()), str(cell_types_path), "cell_type")

        payload: dict[str, object] = {
            "mixture_path": str(mixture_path),
            "scrna_path": str(scrna_path),
            "cell_types_path": str(cell_types_path),
            "tum_key": tum_key,
            "update_gibbs": update_gibbs,
            "n_cores": n_cores,
            "which_theta": which_theta,
            "state_or_type": state_or_type,
            "outlier_cut": outlier_cut,
            "outlier_fraction": outlier_fraction,
            "pseudo_min": pseudo_min,
        }

        if cell_states is not None:
            cell_states_path = workdir_path / "cell_states.tsv"
            write_vector(pd.Series(cell_states.to_numpy()), str(cell_states_path), "cell_state")
            payload["cell_states_path"] = str(cell_states_path)

        result, _ = run_r_tool(tool="bayesprism", payload=payload)

    return DeconvolutionResult(
        tool="bayesprism",
        proportions=result,
        metadata={
            "mixture_transform": mixture_transform,
            "scrna_transform": scrna_transform,
            "update_gibbs": update_gibbs,
            "which_theta": which_theta,
            "state_or_type": state_or_type,
            "outlier_cut": outlier_cut,
            "outlier_fraction": outlier_fraction,
            "pseudo_min": pseudo_min,
        },
    )
