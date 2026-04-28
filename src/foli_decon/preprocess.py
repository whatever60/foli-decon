"""Shared preprocessing utilities across deconvolution wrappers."""

from collections.abc import Sequence

import numpy as np
import pandas as pd

VALID_TRANSFORMS = ("raw", "cpm", "cpm_log1p", "tpm", "tpm_log1p")


def validate_expression_matrix(matrix: pd.DataFrame, matrix_name: str) -> None:
    """Validate that a matrix is non-empty, numeric, and uniquely indexed."""

    if matrix.empty:
        raise ValueError(f"{matrix_name} is empty")
    if not matrix.index.is_unique:
        raise ValueError(f"{matrix_name} has duplicated gene names")
    if not matrix.columns.is_unique:
        raise ValueError(f"{matrix_name} has duplicated sample names")
    is_numeric = matrix.dtypes.apply(lambda dtype: np.issubdtype(dtype, np.number)).all()
    if not is_numeric:
        raise ValueError(f"{matrix_name} contains non-numeric values")


def normalize_gene_names(matrix: pd.DataFrame) -> pd.DataFrame:
    """Standardize gene names by replacing dots with dashes."""

    normalized = matrix.copy()
    normalized.index = normalized.index.astype(str).str.replace(".", "-", regex=False)
    return normalized


def intersect_gene_sets(matrices: Sequence[pd.DataFrame]) -> list[pd.DataFrame]:
    """Subset matrices to their shared gene intersection in original order."""

    shared_genes = set(matrices[0].index)
    for matrix in matrices[1:]:
        shared_genes = shared_genes.intersection(matrix.index)
    ordered_shared_genes = [gene for gene in matrices[0].index if gene in shared_genes]
    intersected = [matrix.loc[ordered_shared_genes] for matrix in matrices]
    return intersected


def counts_to_cpm(counts: pd.DataFrame) -> pd.DataFrame:
    """Convert raw counts to counts per million (CPM)."""

    validate_expression_matrix(counts, "counts")
    library_sizes = counts.sum(axis=0)
    if (library_sizes <= 0).any():
        raise ValueError("All library sizes must be positive for CPM conversion")
    return counts.div(library_sizes, axis=1) * 1_000_000.0


def counts_to_tpm(counts: pd.DataFrame, gene_lengths: pd.Series) -> pd.DataFrame:
    """Convert raw counts to transcripts per million (TPM)."""

    validate_expression_matrix(counts, "counts")
    gene_lengths = gene_lengths.loc[counts.index]
    if (gene_lengths <= 0).any():
        raise ValueError("All gene lengths must be positive for TPM conversion")
    lengths_kb = gene_lengths / 1_000.0
    rpk = counts.div(lengths_kb, axis=0)
    scaling = rpk.sum(axis=0) / 1_000_000.0
    if (scaling <= 0).any():
        raise ValueError("All TPM scaling values must be positive")
    return rpk.div(scaling, axis=1)


def log1p_transform(matrix: pd.DataFrame) -> pd.DataFrame:
    """Apply element-wise natural log(1+x) transformation."""

    return np.log1p(matrix)


def apply_transform(
    matrix: pd.DataFrame,
    transform: str,
    gene_lengths: pd.Series | None = None,
) -> pd.DataFrame:
    """Apply a named transformation to an expression matrix."""

    if transform not in VALID_TRANSFORMS:
        raise ValueError(f"transform must be one of {VALID_TRANSFORMS}")

    if transform == "raw":
        return matrix.copy()
    if transform == "cpm":
        return counts_to_cpm(matrix)
    if transform == "cpm_log1p":
        return log1p_transform(counts_to_cpm(matrix))
    if transform == "tpm":
        if gene_lengths is None:
            raise ValueError("gene_lengths are required for TPM transformation")
        return counts_to_tpm(matrix, gene_lengths)
    if transform == "tpm_log1p":
        if gene_lengths is None:
            raise ValueError("gene_lengths are required for TPM transformation")
        return log1p_transform(counts_to_tpm(matrix, gene_lengths))

    raise ValueError(f"Unsupported transform: {transform}")


def build_signature_from_scrna(
    scrna_counts: pd.DataFrame,
    cell_types: pd.Series,
    transform: str = "cpm",
    gene_lengths: pd.Series | None = None,
    min_cells_per_type: int = 1,
) -> pd.DataFrame:
    """Aggregate single-cell expression to a gene-by-cell-type signature matrix."""

    validate_expression_matrix(scrna_counts, "scrna_counts")
    if len(cell_types) != scrna_counts.shape[1]:
        raise ValueError("cell_types length must match number of single-cell columns")

    cell_types_aligned = pd.Series(cell_types.to_numpy(), index=scrna_counts.columns)
    cell_type_counts = cell_types_aligned.value_counts()
    if (cell_type_counts < min_cells_per_type).any():
        raise ValueError("At least one cell type has fewer cells than min_cells_per_type")

    transformed = apply_transform(scrna_counts, transform=transform, gene_lengths=gene_lengths)
    signature = transformed.T.groupby(cell_types_aligned).mean().T
    signature.columns.name = "cell_type"
    return signature


def align_mixture_and_signature(
    mixture: pd.DataFrame,
    signature: pd.DataFrame,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Normalize gene names and align mixture/signature matrices by gene intersection."""

    mixture_norm = normalize_gene_names(mixture)
    signature_norm = normalize_gene_names(signature)
    mixture_aligned, signature_aligned = intersect_gene_sets([mixture_norm, signature_norm])
    return mixture_aligned, signature_aligned


def align_mixture_and_scrna(
    mixture: pd.DataFrame,
    scrna_counts: pd.DataFrame,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Normalize gene names and align mixture/single-cell matrices by gene intersection."""

    mixture_norm = normalize_gene_names(mixture)
    scrna_norm = normalize_gene_names(scrna_counts)
    mixture_aligned, scrna_aligned = intersect_gene_sets([mixture_norm, scrna_norm])
    return mixture_aligned, scrna_aligned
