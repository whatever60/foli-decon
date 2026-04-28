"""Unit tests for two-split benchmark workflow."""

import numpy as np
import pandas as pd
from anndata import AnnData

from foli_decon.benchmark import run_two_split_benchmark_from_adata
from foli_decon.models import DeconvolutionResult


def _build_tiny_adata() -> AnnData:
    """Create a small synthetic AnnData object with raw_counts layer."""

    rng = np.random.default_rng(7)
    counts = rng.poisson(lam=5, size=(24, 12)).astype(float)
    obs = pd.DataFrame(
        {
            "celltype": ["A"] * 8 + ["B"] * 8 + ["C"] * 8,
            "donor": ["d1", "d2"] * 12,
        },
        index=[f"cell_{i}" for i in range(24)],
    )
    var = pd.DataFrame(index=[f"Gene_{i}" for i in range(12)])
    adata = AnnData(X=counts, obs=obs, var=var)
    adata.layers["raw_counts"] = counts
    return adata


def test_two_split_benchmark_returns_metrics(monkeypatch) -> None:
    """Two-split benchmark should emit one row per tool and scenario."""

    def fake_run_deconvolution(tool: str, **kwargs) -> DeconvolutionResult:
        """Return simple uniform predictions aligned to provided references."""

        mixture = kwargs["mixture"]
        if "signature" in kwargs:
            output_columns = kwargs["signature"].columns.tolist()
        elif "cell_types" in kwargs:
            output_columns = pd.Index(kwargs["cell_types"]).unique().tolist()
        else:
            output_columns = ["A", "B", "C"]
        values = np.ones((mixture.shape[1], len(output_columns)), dtype=float)
        values = values / values.sum(axis=1, keepdims=True)
        prediction = pd.DataFrame(values, index=mixture.columns, columns=output_columns)
        return DeconvolutionResult(tool=tool, proportions=prediction, metadata={})

    monkeypatch.setattr("foli_decon.benchmark.run_deconvolution", fake_run_deconvolution)

    adata = _build_tiny_adata()
    result = run_two_split_benchmark_from_adata(
        adata=adata,
        dataset_name="tiny",
        tools=["deconrnaseq", "music"],
        cell_type_col="celltype",
        batch_col="donor",
        count_layer="raw_counts",
        min_cells_per_type=4,
        max_reference_cells_per_type=4,
        max_evaluation_cells_per_type=4,
        n_samples_per_scenario=3,
        cells_per_bulk=30,
        seed=11,
    )

    assert result.metrics.shape[0] == 6
    assert set(result.metrics["scenario"].tolist()) == {
        "dirichlet_uniform",
        "one_type_dominant",
        "one_cell_type_dominate",
    }
    assert set(result.metrics["tool"].tolist()) == {"deconrnaseq", "music"}
    assert (result.metrics["status"] == "ok").all()
    assert (result.metrics["n_shared_cell_types"] > 0).all()


def test_two_split_benchmark_handles_gene_name_collisions(monkeypatch) -> None:
    """Benchmark preprocessing should handle dot-to-dash gene collisions."""

    def fake_run_deconvolution(tool: str, **kwargs) -> DeconvolutionResult:
        """Return simple predictions for all mixture samples."""

        mixture = kwargs["mixture"]
        output_columns = kwargs["signature"].columns.tolist()
        values = np.ones((mixture.shape[1], len(output_columns)), dtype=float)
        values = values / values.sum(axis=1, keepdims=True)
        prediction = pd.DataFrame(values, index=mixture.columns, columns=output_columns)
        return DeconvolutionResult(tool=tool, proportions=prediction, metadata={})

    monkeypatch.setattr("foli_decon.benchmark.run_deconvolution", fake_run_deconvolution)

    adata = _build_tiny_adata()
    adata.var_names = pd.Index(
        [
            "Gene.0",
            "Gene-0",
            "Gene.1",
            "Gene.2",
            "Gene.3",
            "Gene.4",
            "Gene.5",
            "Gene.6",
            "Gene.7",
            "Gene.8",
            "Gene.9",
            "Gene.10",
        ]
    )
    result = run_two_split_benchmark_from_adata(
        adata=adata,
        dataset_name="tiny",
        tools=["deconrnaseq"],
        cell_type_col="celltype",
        batch_col="donor",
        count_layer="raw_counts",
        min_cells_per_type=4,
        max_reference_cells_per_type=4,
        max_evaluation_cells_per_type=4,
        n_samples_per_scenario=2,
        cells_per_bulk=20,
        seed=13,
    )

    assert (result.metrics["status"] == "ok").all()


def test_two_split_benchmark_key_stratified_split(monkeypatch) -> None:
    """Key-stratified split should keep one donor group in reference cells."""

    seen_batch_sets: list[set[str]] = []

    def fake_run_deconvolution(tool: str, **kwargs) -> DeconvolutionResult:
        """Return uniform predictions and capture batch IDs for reference cells."""

        if "batch_ids" in kwargs:
            seen_batch_sets.append(set(kwargs["batch_ids"].astype(str).tolist()))
        mixture = kwargs["mixture"]
        output_columns = pd.Index(kwargs["cell_types"]).unique().tolist()
        values = np.ones((mixture.shape[1], len(output_columns)), dtype=float)
        values = values / values.sum(axis=1, keepdims=True)
        prediction = pd.DataFrame(values, index=mixture.columns, columns=output_columns)
        return DeconvolutionResult(tool=tool, proportions=prediction, metadata={})

    monkeypatch.setattr("foli_decon.benchmark.run_deconvolution", fake_run_deconvolution)

    adata = _build_tiny_adata()
    result = run_two_split_benchmark_from_adata(
        adata=adata,
        dataset_name="tiny",
        tools=["music"],
        cell_type_col="celltype",
        batch_col="donor",
        count_layer="raw_counts",
        min_cells_per_type=4,
        max_reference_cells_per_type=4,
        max_evaluation_cells_per_type=4,
        n_samples_per_scenario=2,
        cells_per_bulk=20,
        split_strategy="key_stratified",
        split_key_col="donor",
        seed=19,
    )

    assert result.metrics.shape[0] == 3
    assert (result.metrics["status"] == "ok").all()
    assert len(seen_batch_sets) == 3
    assert all(len(batch_set) == 1 for batch_set in seen_batch_sets)
