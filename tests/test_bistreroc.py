"""Focused tests for the BistreRoc sidecar adapter."""

import gzip
from pathlib import Path

import pandas as pd
import pytest

from foli_decon import (
    run_bistreroc,
    run_deconvolution,
    supported_reference_trainers,
    supported_tools,
    train_bistreroc_reference,
    train_reference,
)
from foli_decon.models import DeconvolutionResult, ReferenceTrainingResult


def test_run_bistreroc_stages_compressed_inputs_and_parses_outputs(
    monkeypatch,
    tmp_path: Path,
) -> None:
    """The adapter should align inputs and separate fractions from diagnostics."""

    captured: dict[str, object] = {}

    def fake_run_command(
        command: list[str],
        cwd: str | None = None,
        log_path: str | None = None,
    ) -> None:
        """Capture staged inputs and write a deterministic native result."""

        mixture_path = Path(command[command.index("--mixture") + 1])
        signature_path = Path(command[command.index("--signature") + 1])
        native_output_path = Path(command[command.index("--output") + 1])
        captured["command"] = command
        captured["mixture"] = pd.read_csv(mixture_path, sep="\t", index_col=0)
        captured["signature"] = pd.read_csv(signature_path, sep="\t", index_col=0)
        pd.DataFrame(
            {
                "T_cell": [0.75, 0.20],
                "B_cell": [0.25, 0.80],
                "P-value": [9999.0, 9999.0],
                "Correlation": [0.96, 0.91],
                "RMSE": [0.04, 0.09],
            },
            index=["sample_1", "sample_2"],
        ).to_csv(native_output_path, sep="\t", index_label="Mixture")

    monkeypatch.setattr(
        "foli_decon.tools.bistreroc.run_command",
        fake_run_command,
    )
    mixture = pd.DataFrame(
        {
            "sample_1": [10.0, 20.0, 30.0],
            "sample_2": [12.0, 22.0, 32.0],
        },
        index=["Gene.A", "GeneB", "mixture_only"],
    )
    signature = pd.DataFrame(
        {
            "T_cell": [50.0, 10.0, 5.0],
            "B_cell": [5.0, 45.0, 10.0],
        },
        index=["Gene-A", "GeneB", "signature_only"],
    )
    output_path = tmp_path / "results" / "proportion.tsv.gz"

    result = run_bistreroc(
        mixture=mixture,
        signature=signature,
        cores=4,
        timeout_seconds=90,
        output_path=str(output_path),
    )

    command = captured["command"]
    assert command[:5] == [
        "mamba",
        "run",
        "-n",
        "foli-decon-bistreroc",
        "bistreroc",
    ]
    assert command[-4:] == ["--cores", "4", "--timeout-seconds", "90"]
    assert captured["mixture"].index.tolist() == ["Gene-A", "GeneB"]
    assert captured["signature"].index.tolist() == ["Gene-A", "GeneB"]
    assert result.tool == "bistreroc"
    assert result.proportion.columns.tolist() == ["T_cell", "B_cell"]
    assert result.score.columns.tolist() == ["Correlation", "RMSE"]
    assert result.p_value.columns.tolist() == ["P-value"]
    assert result.metadata["shared_gene_count"] == 2
    assert result.output_paths["proportion"] == str(output_path.resolve())
    pd.testing.assert_frame_equal(
        pd.read_csv(output_path, sep="\t", index_col=0),
        result.proportion,
    )


def test_train_bistreroc_reference_stages_integer_counts_with_cell_type_labels(
    monkeypatch,
    tmp_path: Path,
) -> None:
    """Reference training should preserve raw counts and use labels as headers."""

    captured: dict[str, object] = {}

    def fake_run_command(
        command: list[str],
        cwd: str | None = None,
        log_path: str | None = None,
    ) -> None:
        """Capture the labeled reference and write a small signature matrix."""

        reference_path = Path(command[command.index("--reference") + 1])
        output_path = Path(command[command.index("--output") + 1])
        with gzip.open(reference_path, "rt", encoding="utf-8") as input_handle:
            captured["reference_text"] = input_handle.read()
        captured["command"] = command
        pd.DataFrame(
            {"T_cell": [100.0, 10.0], "B_cell": [5.0, 90.0]},
            index=["Gene-1", "Gene2"],
        ).to_csv(output_path, sep="\t", index_label="Gene")

    monkeypatch.setattr(
        "foli_decon.tools.bistreroc.run_command",
        fake_run_command,
    )
    counts = pd.DataFrame(
        {
            "cell_1": [1.0, 3.0],
            "cell_2": [2.0, 4.0],
            "cell_3": [5.0, 6.0],
        },
        index=["Gene.1", "Gene2"],
    )
    cell_types = pd.Series(["T_cell", "T_cell", "B_cell"])
    output_path = tmp_path / "signature.tsv"

    result = train_bistreroc_reference(
        scrna_counts=counts,
        cell_types=cell_types,
        output_path=str(output_path),
        minimum_markers_per_cell_type=7,
        maximum_markers_per_cell_type=11,
        q_value_threshold=0.05,
        cores=3,
    )

    command = captured["command"]
    assert command[5] == "build-signature"
    assert command[command.index("--minimum-markers-per-cell-type") + 1] == "7"
    assert command[command.index("--maximum-markers-per-cell-type") + 1] == "11"
    assert command[command.index("--q-value-threshold") + 1] == "0.05"
    reference_lines = captured["reference_text"].splitlines()
    assert reference_lines[0].split("\t")[1:] == ["T_cell", "T_cell", "B_cell"]
    assert reference_lines[1] == "Gene-1\t1\t2\t5"
    assert result.tool == "bistreroc"
    assert result.signature_path == str(output_path.resolve())
    assert result.output_paths["signature"] == str(output_path.resolve())
    assert result.signature.columns.tolist() == ["T_cell", "B_cell"]
    assert result.metadata["n_cell_types"] == 2


def test_train_bistreroc_reference_rejects_noninteger_counts(tmp_path: Path) -> None:
    """BistreRoc signature generation should receive raw integer counts."""

    counts = pd.DataFrame(
        {"cell_1": [1.5, 2.0]},
        index=["Gene1", "Gene2"],
    )
    with pytest.raises(ValueError, match="requires integer raw counts"):
        train_bistreroc_reference(
            scrna_counts=counts,
            cell_types=pd.Series(["T_cell"]),
            output_path=str(tmp_path / "signature.tsv"),
        )


def test_bistreroc_is_registered_in_both_public_dispatchers(monkeypatch) -> None:
    """The unified APIs should expose both BistreRoc computational stages."""

    deconvolution_result = DeconvolutionResult(tool="bistreroc")
    reference_result = ReferenceTrainingResult(tool="bistreroc")
    monkeypatch.setattr(
        "foli_decon.api.run_bistreroc",
        lambda **kwargs: deconvolution_result,
    )
    monkeypatch.setattr(
        "foli_decon.api.train_bistreroc_reference",
        lambda **kwargs: reference_result,
    )

    assert "bistreroc" in supported_tools()
    assert "bistreroc" in supported_reference_trainers()
    assert run_deconvolution("BistreRoc", sentinel=True) is deconvolution_result
    assert train_reference("bistreroc", sentinel=True) is reference_result


def test_bistreroc_sidecar_is_pinned_to_the_release_tag() -> None:
    """The temporary integration should install the immutable release tag."""

    environment_text = Path("environment-bistreroc.yml").read_text(encoding="utf-8")
    assert "https://github.com/whatever60/BistreRoc.git@v0.1.0" in environment_text
