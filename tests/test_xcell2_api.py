"""Unit tests for xCell2 training and end-to-end scoring APIs."""

from pathlib import Path

import pandas as pd

from foli_decon import run_deconvolution, train_reference


def test_train_xcell2_reference_payload(monkeypatch, tmp_path) -> None:
    """xCell2 reference training should write a reference object through the R wrapper."""

    calls = []

    def fake_run_r_tool(tool: str, payload: dict[str, object], **kwargs):
        """Capture training payload and create a fake reference object."""

        calls.append((tool, payload, kwargs))
        Path(str(payload["xcell2_object_path"])).write_text("fake reference", encoding="utf-8")
        return pd.DataFrame({"status": ["ok"]}, index=["xcell2_train"]), payload

    monkeypatch.setattr("foli_decon.tools.bulk_tools.run_r_tool", fake_run_r_tool)

    scrna_counts = pd.DataFrame(
        {
            "cell_1": [10.0, 0.0],
            "cell_2": [30.0, 10.0],
            "cell_3": [0.0, 20.0],
            "cell_4": [5.0, 25.0],
        },
        index=["Gene.1", "Gene2"],
    )
    cell_types = pd.Series(["A", "A", "B", "B"], index=scrna_counts.columns)
    output_path = tmp_path / "xcell2_reference.rds"

    result = train_reference(
        tool="xcell2",
        scrna_counts=scrna_counts,
        cell_types=cell_types,
        output_path=str(output_path),
        return_signatures=True,
        min_sc_genes=2,
        min_pb_cells=1,
        min_pb_samples=1,
        xcell2_workers=2,
        log_path=str(tmp_path / "train.log"),
    )

    assert result.rds_path == str(output_path)
    assert result.metadata["reference_mode"] == "signature_only"
    assert result.metadata["raw_scores"] is True
    assert result.metadata["n_cell_types"] == 2
    assert calls[0][0] == "xcell2_train"
    assert calls[0][1]["return_signatures"] is True
    assert calls[0][1]["xcell2_workers"] == 2


def test_run_xcell2_can_train_then_score(monkeypatch, tmp_path) -> None:
    """run_deconvolution(tool='xcell2') should train from scRNA-seq when no object path is supplied."""

    calls = []

    def fake_run_r_tool(tool: str, payload: dict[str, object], **kwargs):
        """Return fake training and scoring outputs."""

        calls.append((tool, payload, kwargs))
        if tool == "xcell2_train":
            Path(str(payload["xcell2_object_path"])).write_text("fake reference", encoding="utf-8")
            return pd.DataFrame({"status": ["ok"]}, index=["xcell2_train"]), payload
        scores = pd.DataFrame({"A": [0.8], "B": [0.2]}, index=["bulk_1"])
        return scores, payload

    monkeypatch.setattr("foli_decon.tools.bulk_tools.run_r_tool", fake_run_r_tool)

    mixture = pd.DataFrame({"bulk_1": [100.0, 50.0]}, index=["Gene1", "Gene2"])
    scrna_counts = pd.DataFrame(
        {
            "cell_1": [10.0, 0.0],
            "cell_2": [30.0, 10.0],
            "cell_3": [0.0, 20.0],
            "cell_4": [5.0, 25.0],
        },
        index=["Gene1", "Gene2"],
    )
    cell_types = pd.Series(["A", "A", "B", "B"], index=scrna_counts.columns)
    reference_path = tmp_path / "persisted_reference.rds"

    result = run_deconvolution(
        tool="xcell2",
        mixture=mixture,
        scrna_counts=scrna_counts,
        cell_types=cell_types,
        xcell2_reference_output_path=str(reference_path),
        return_signatures=True,
        transform="raw",
        min_shared_genes=0.1,
        min_sc_genes=2,
        min_pb_cells=1,
        min_pb_samples=1,
    )

    assert result.score.loc["bulk_1", "A"] == 0.8
    assert result.metadata["raw_scores"] is True
    assert result.metadata["spillover"] is False
    assert result.metadata["training_reference_mode"] == "signature_only"
    assert result.metadata["reference_persisted"] is True
    assert [call[0] for call in calls] == ["xcell2_train", "xcell2"]
