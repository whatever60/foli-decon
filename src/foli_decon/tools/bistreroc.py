"""Adapters for the BistreRoc command-line package."""

import tempfile
from pathlib import Path

import numpy as np
import pandas as pd

from foli_decon.io import (
    read_expression_matrix,
    read_sample_by_celltype,
    write_expression_matrix,
)
from foli_decon.models import DeconvolutionResult, ReferenceTrainingResult
from foli_decon.preprocess import (
    apply_transform,
    intersect_gene_sets,
    normalize_gene_names,
    validate_expression_matrix,
)
from foli_decon.python_runner import run_command


def _bistreroc_command_prefix(bistreroc_env: str | None) -> list[str]:
    """Build the command prefix for an installed or sidecar BistreRoc CLI."""

    if bistreroc_env is None:
        return ["bistreroc"]
    return ["mamba", "run", "-n", bistreroc_env, "bistreroc"]


def run_bistreroc(
    mixture: pd.DataFrame,
    signature: pd.DataFrame,
    mixture_transform: str = "raw",
    signature_transform: str = "raw",
    gene_lengths: pd.Series | None = None,
    bistreroc_env: str | None = "foli-decon-bistreroc",
    cores: int = 1,
    timeout_seconds: int = 3600,
    output_path: str | None = None,
    log_path: str | None = None,
) -> DeconvolutionResult:
    """Estimate cell fractions with BistreRoc in its dedicated sidecar."""

    validate_expression_matrix(mixture, "mixture")
    validate_expression_matrix(signature, "signature")
    mixture_transformed = apply_transform(
        normalize_gene_names(mixture),
        transform=mixture_transform,
        gene_lengths=gene_lengths,
    )
    signature_transformed = apply_transform(
        normalize_gene_names(signature),
        transform=signature_transform,
        gene_lengths=gene_lengths,
    )
    mixture_aligned, signature_aligned = intersect_gene_sets(
        [mixture_transformed, signature_transformed]
    )

    with tempfile.TemporaryDirectory(prefix="foli_bistreroc_") as workdir:
        workdir_path = Path(workdir)
        mixture_path = workdir_path / "mixture.tsv.gz"
        signature_path = workdir_path / "signature.tsv.gz"
        native_output_path = workdir_path / "fractions.tsv"
        write_expression_matrix(mixture_aligned, str(mixture_path))
        write_expression_matrix(signature_aligned, str(signature_path))
        command = _bistreroc_command_prefix(bistreroc_env) + [
            "deconvolve",
            "--mixture",
            str(mixture_path),
            "--signature",
            str(signature_path),
            "--output",
            str(native_output_path),
            "--cores",
            str(cores),
            "--timeout-seconds",
            str(timeout_seconds),
        ]
        run_command(command, log_path=log_path)
        native_result = read_sample_by_celltype(str(native_output_path))

    p_value = native_result[["P-value"]]
    score = native_result[["Correlation", "RMSE"]]
    proportion = native_result.drop(columns=["P-value", "Correlation", "RMSE"])
    output_paths: dict[str, str] = {}
    if output_path is not None:
        resolved_output_path = Path(output_path).resolve()
        resolved_output_path.parent.mkdir(parents=True, exist_ok=True)
        proportion.to_csv(resolved_output_path, sep="\t")
        output_paths["proportion"] = str(resolved_output_path)
    if log_path is not None:
        output_paths["log"] = str(Path(log_path).resolve())

    return DeconvolutionResult(
        tool="bistreroc",
        proportion=proportion,
        score=score,
        p_value=p_value,
        metadata={
            "mixture_transform": mixture_transform,
            "signature_transform": signature_transform,
            "bistreroc_env": "" if bistreroc_env is None else bistreroc_env,
            "cores": cores,
            "timeout_seconds": timeout_seconds,
            "mixture_input_gene_count": mixture.shape[0],
            "signature_input_gene_count": signature.shape[0],
            "shared_gene_count": mixture_aligned.shape[0],
            "sample_count": mixture_aligned.shape[1],
            "cell_type_count": signature_aligned.shape[1],
        },
        output_paths=output_paths,
    )


def train_bistreroc_reference(
    scrna_counts: pd.DataFrame,
    cell_types: pd.Series,
    output_path: str,
    minimum_markers_per_cell_type: int = 300,
    maximum_markers_per_cell_type: int = 500,
    q_value_threshold: float = 0.01,
    bistreroc_env: str | None = "foli-decon-bistreroc",
    cores: int = 1,
    timeout_seconds: int = 3600,
    log_path: str | None = None,
) -> ReferenceTrainingResult:
    """Build and persist a BistreRoc signature from labeled raw counts."""

    validate_expression_matrix(scrna_counts, "scrna_counts")
    if len(cell_types) != scrna_counts.shape[1]:
        raise ValueError("cell_types length must match number of single-cell columns")
    normalized_counts = normalize_gene_names(scrna_counts)
    count_values = normalized_counts.to_numpy(dtype=np.float64)
    if not np.isfinite(count_values).all():
        raise ValueError("scrna_counts must contain only finite values")
    if (count_values < 0).any():
        raise ValueError("scrna_counts must be nonnegative")
    if not np.equal(count_values, np.floor(count_values)).all():
        raise ValueError("BistreRoc signature generation requires integer raw counts")
    cell_type_labels = pd.Series(cell_types.to_numpy())
    if cell_type_labels.isna().any():
        raise ValueError("cell_types cannot contain missing values")
    staged_reference = normalized_counts.astype(np.int64)
    staged_reference.columns = cell_type_labels.astype(str).to_numpy()

    resolved_output_path = Path(output_path).resolve()
    resolved_output_path.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(prefix="foli_bistreroc_signature_") as workdir:
        reference_path = Path(workdir) / "reference.tsv.gz"
        write_expression_matrix(staged_reference, str(reference_path))
        command = _bistreroc_command_prefix(bistreroc_env) + [
            "build-signature",
            "--reference",
            str(reference_path),
            "--output",
            str(resolved_output_path),
            "--minimum-markers-per-cell-type",
            str(minimum_markers_per_cell_type),
            "--maximum-markers-per-cell-type",
            str(maximum_markers_per_cell_type),
            "--q-value-threshold",
            str(q_value_threshold),
            "--cores",
            str(cores),
            "--timeout-seconds",
            str(timeout_seconds),
        ]
        run_command(command, log_path=log_path)

    signature = read_expression_matrix(str(resolved_output_path))
    output_paths = {"signature": str(resolved_output_path)}
    if log_path is not None:
        output_paths["log"] = str(Path(log_path).resolve())
    return ReferenceTrainingResult(
        tool="bistreroc",
        signature=signature,
        signature_path=str(resolved_output_path),
        metadata={
            "bistreroc_env": "" if bistreroc_env is None else bistreroc_env,
            "minimum_markers_per_cell_type": minimum_markers_per_cell_type,
            "maximum_markers_per_cell_type": maximum_markers_per_cell_type,
            "q_value_threshold": q_value_threshold,
            "cores": cores,
            "timeout_seconds": timeout_seconds,
            "n_cells": staged_reference.shape[1],
            "n_genes": staged_reference.shape[0],
            "n_cell_types": cell_type_labels.nunique(),
        },
        output_paths=output_paths,
    )
