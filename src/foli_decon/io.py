"""Lightweight file I/O helpers for expression data."""

from pathlib import Path

import pandas as pd


def read_expression_matrix(path: str, sep: str = "\t") -> pd.DataFrame:
    """Read an expression matrix with genes as index and samples as columns."""

    matrix = pd.read_csv(path, sep=sep, index_col=0)
    matrix.index = matrix.index.astype(str)
    matrix.columns = matrix.columns.astype(str)
    return matrix


def write_expression_matrix(matrix: pd.DataFrame, path: str, sep: str = "\t") -> None:
    """Write an expression matrix with row names preserved."""

    output_path = Path(path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    matrix.to_csv(output_path, sep=sep)


def write_vector(vector: pd.Series, path: str, column_name: str) -> None:
    """Write a one-column vector table for R-side metadata inputs."""

    output_path = Path(path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    vector.to_frame(name=column_name).to_csv(output_path, sep="\t", index=False)


def read_sample_by_celltype(path: str) -> pd.DataFrame:
    """Read wrapper outputs formatted as sample-by-celltype tables."""

    result = pd.read_csv(path, sep="\t", index_col=0)
    result.index = result.index.astype(str)
    result.columns = result.columns.astype(str)
    return result
