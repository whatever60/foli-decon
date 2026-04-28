"""Unit tests for preprocessing utilities."""

import numpy as np
import pandas as pd

from foli_decon.preprocess import (
    align_mixture_and_signature,
    build_signature_from_scrna,
    counts_to_cpm,
)


def test_counts_to_cpm_column_sums_are_one_million() -> None:
    """CPM conversion should normalize each sample column to 1e6."""

    counts = pd.DataFrame(
        {
            "sample_a": [10.0, 30.0],
            "sample_b": [5.0, 5.0],
        },
        index=["Gene1", "Gene2"],
    )
    cpm = counts_to_cpm(counts)
    np.testing.assert_allclose(cpm.sum(axis=0).to_numpy(), np.array([1_000_000.0, 1_000_000.0]))


def test_build_signature_from_scrna_groups_cells() -> None:
    """Signature generation should average transformed cells by cell type."""

    scrna = pd.DataFrame(
        {
            "cell_1": [10.0, 0.0],
            "cell_2": [30.0, 10.0],
            "cell_3": [20.0, 20.0],
        },
        index=["Gene1", "Gene2"],
    )
    cell_types = pd.Series(["A", "A", "B"])
    signature = build_signature_from_scrna(scrna, cell_types, transform="raw")
    assert list(signature.columns) == ["A", "B"]
    np.testing.assert_allclose(signature.loc["Gene1", "A"], 20.0)
    np.testing.assert_allclose(signature.loc["Gene1", "B"], 20.0)


def test_align_mixture_and_signature_intersects_genes() -> None:
    """Alignment should keep only shared genes in mixture order."""

    mixture = pd.DataFrame({"s1": [1.0, 2.0, 3.0]}, index=["A", "B", "C"])
    signature = pd.DataFrame({"ct1": [10.0, 20.0]}, index=["C", "A"])
    mixture_aligned, signature_aligned = align_mixture_and_signature(mixture, signature)
    assert list(mixture_aligned.index) == ["A", "C"]
    assert list(signature_aligned.index) == ["A", "C"]
