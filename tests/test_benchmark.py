"""Unit tests for two-split benchmark workflow."""

import runpy

import numpy as np
import pandas as pd
from anndata import AnnData

from foli_decon.benchmark import run_two_split_benchmark_from_adata
from foli_decon.models import DeconvolutionResult, ReferenceTrainingResult


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
        return DeconvolutionResult(tool=tool, proportion=prediction, metadata={})

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


def test_cli_default_excludes_manual_or_sidecar_methods() -> None:
    """CLI defaults should skip tools needing credentials, marker sets, or sidecar runtimes."""

    script_globals = runpy.run_path("scripts/run_two_split_benchmark.py")
    default_excluded_tools = script_globals["BENCHMARK_DEFAULT_EXCLUDED_TOOLS"]
    expected_exclusions = {
        "autogenes",
        "bistreroc",
        "blade",
        "blue",
        "cdseq",
        "cibersortx",
        "dissect",
        "instaprism",
        "music2",
        "psea",
        "scaden",
        "scdc",
        "tape",
    }

    assert expected_exclusions.issubset(default_excluded_tools)


def test_two_split_benchmark_passes_references_to_explicit_sidecar_tools(monkeypatch) -> None:
    """Explicit sidecar tools should receive the reference input type their wrappers need."""

    seen_keys: dict[str, set[str]] = {}

    def fake_run_deconvolution(tool: str, **kwargs) -> DeconvolutionResult:
        """Capture reference kwargs and return uniform predictions."""

        seen_keys[tool] = set(kwargs)
        mixture = kwargs["mixture"]
        if tool in {"autogenes", "bistreroc"}:
            output_columns = kwargs["signature"].columns.tolist()
        else:
            output_columns = pd.Index(kwargs["cell_types"]).unique().tolist()
            assert "scrna_counts" in kwargs
        values = np.ones((mixture.shape[1], len(output_columns)), dtype=float)
        values = values / values.sum(axis=1, keepdims=True)
        prediction = pd.DataFrame(values, index=mixture.columns, columns=output_columns)
        return DeconvolutionResult(tool=tool, proportion=prediction, metadata={})

    monkeypatch.setattr("foli_decon.benchmark.run_deconvolution", fake_run_deconvolution)

    adata = _build_tiny_adata()
    result = run_two_split_benchmark_from_adata(
        adata=adata,
        dataset_name="tiny",
        tools=["autogenes", "bistreroc", "blade", "blue", "tape"],
        cell_type_col="celltype",
        batch_col="donor",
        count_layer="raw_counts",
        min_cells_per_type=4,
        max_reference_cells_per_type=4,
        max_evaluation_cells_per_type=4,
        n_samples_per_scenario=2,
        cells_per_bulk=20,
        seed=12,
    )

    assert result.metrics.shape[0] == 15
    assert (result.metrics["status"] == "ok").all()
    assert "signature" in seen_keys["autogenes"]
    assert "signature" in seen_keys["bistreroc"]
    assert {"scrna_counts", "cell_types"}.issubset(seen_keys["blade"])
    assert {"scrna_counts", "cell_types", "batch_ids"}.issubset(seen_keys["blue"])
    assert {"scrna_counts", "cell_types"}.issubset(seen_keys["tape"])


def test_two_split_benchmark_handles_gene_name_collisions(monkeypatch) -> None:
    """Benchmark preprocessing should handle dot-to-dash gene collisions."""

    def fake_run_deconvolution(tool: str, **kwargs) -> DeconvolutionResult:
        """Return simple predictions for all mixture samples."""

        mixture = kwargs["mixture"]
        output_columns = kwargs["signature"].columns.tolist()
        values = np.ones((mixture.shape[1], len(output_columns)), dtype=float)
        values = values / values.sum(axis=1, keepdims=True)
        prediction = pd.DataFrame(values, index=mixture.columns, columns=output_columns)
        return DeconvolutionResult(tool=tool, proportion=prediction, metadata={})

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
        return DeconvolutionResult(tool=tool, proportion=prediction, metadata={})

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


def test_two_split_benchmark_trains_xcell2_reference_for_scores(monkeypatch) -> None:
    """xCell2 benchmark should train a split-specific reference and report score metrics."""

    trained_paths: list[str] = []

    def fake_train_reference(tool: str, **kwargs) -> ReferenceTrainingResult:
        """Capture xCell2 reference training calls."""

        output_path = kwargs["output_path"]
        trained_paths.append(output_path)
        return ReferenceTrainingResult(
            tool=tool,
            rds_path=output_path,
            metadata={"return_signatures": kwargs["return_signatures"]},
        )

    def fake_run_deconvolution(tool: str, **kwargs) -> DeconvolutionResult:
        """Return score-like xCell2 outputs aligned to the tiny cell types."""

        mixture = kwargs["mixture"]
        assert kwargs["xcell2_object_path"] == trained_paths[0]
        assert kwargs["raw_scores"] is True
        assert kwargs["spillover"] is False
        columns = ["A", "B", "C"]
        sample_values = np.arange(1, mixture.shape[1] + 1, dtype=float)[:, None]
        type_values = np.arange(1, len(columns) + 1, dtype=float)[None, :]
        prediction = pd.DataFrame(sample_values * type_values, index=mixture.columns, columns=columns)
        return DeconvolutionResult(tool=tool, score=prediction, metadata={})

    monkeypatch.setattr("foli_decon.benchmark.train_reference", fake_train_reference)
    monkeypatch.setattr("foli_decon.benchmark.run_deconvolution", fake_run_deconvolution)

    adata = _build_tiny_adata()
    result = run_two_split_benchmark_from_adata(
        adata=adata,
        dataset_name="tiny",
        tools=["xcell2"],
        cell_type_col="celltype",
        batch_col="donor",
        count_layer="raw_counts",
        min_cells_per_type=4,
        max_reference_cells_per_type=4,
        max_evaluation_cells_per_type=4,
        n_samples_per_scenario=3,
        cells_per_bulk=30,
        seed=23,
    )

    assert len(trained_paths) == 1
    assert result.metrics.shape[0] == 3
    assert (result.metrics["status"] == "ok").all()
    assert result.metrics["mae"].isna().all()
    assert result.metrics["score_spearman"].notna().all()
    assert (result.metrics["n_shared_cell_types"] == 3).all()


def test_two_split_benchmark_keeps_score_run_ok_without_label_overlap(monkeypatch) -> None:
    """Score-output tools should stay run-successful when labels cannot be compared directly."""

    def fake_run_deconvolution(tool: str, **kwargs) -> DeconvolutionResult:
        """Return xCell-like scores with labels absent from the reference truth."""

        mixture = kwargs["mixture"]
        prediction = pd.DataFrame(
            {
                "Epithelial cells": [1.0] * mixture.shape[1],
                "Macrophages": [2.0] * mixture.shape[1],
            },
            index=mixture.columns,
        )
        return DeconvolutionResult(tool=tool, score=prediction, metadata={})

    monkeypatch.setattr("foli_decon.benchmark.run_deconvolution", fake_run_deconvolution)

    adata = _build_tiny_adata()
    result = run_two_split_benchmark_from_adata(
        adata=adata,
        dataset_name="tiny",
        tools=["xcell"],
        cell_type_col="celltype",
        batch_col="donor",
        count_layer="raw_counts",
        min_cells_per_type=4,
        max_reference_cells_per_type=4,
        max_evaluation_cells_per_type=4,
        n_samples_per_scenario=2,
        cells_per_bulk=20,
        seed=29,
    )

    assert result.metrics.shape[0] == 3
    assert (result.metrics["status"] == "ok").all()
    assert result.metrics["score_spearman"].isna().all()
    assert result.metrics["error_message"].str.contains("score metrics not computed").all()
