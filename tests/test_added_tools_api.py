"""Unit tests for newly added deconvolution wrapper APIs."""

import json
from pathlib import Path

import pandas as pd
import pytest

from foli_decon import (
    run_deconvolution,
    run_music2,
    select_features,
    select_features_autogenes,
    select_features_markermap,
    select_features_scgenefit,
    supported_feature_selectors,
    supported_tools,
)


def test_run_xcell_rejects_small_targeted_panel(monkeypatch) -> None:
    """xCell should fail before R when fewer than its minimum input genes are provided."""

    def fake_run_r_tool(tool: str, payload: dict[str, object], **kwargs):
        """Fail if the xCell preflight check does not stop execution first."""

        raise AssertionError("xCell R wrapper should not be called for too-small panels")

    monkeypatch.setattr("foli_decon.tools.bulk_tools.run_r_tool", fake_run_r_tool)

    mixture = pd.DataFrame({"bulk_1": [10.0, 30.0]}, index=["Gene1", "Gene2"])
    with pytest.raises(ValueError, match="xCell requires at least 5000 input genes"):
        run_deconvolution(tool="xcell", mixture=mixture, transform="raw")


def test_supported_tools_excludes_deferred_web_or_unresolved_methods() -> None:
    """Deferred or web-only methods should not appear in supported tool names."""

    names = set(supported_tools())
    assert "tape" in names
    assert "blade" in names
    assert "blue" in names
    assert "abis" in names
    assert "estimate" in names
    assert "timer" not in names
    assert "timer2" not in names
    assert "timer3" not in names


def test_run_xcell_large_matrix_payload(monkeypatch) -> None:
    """xCell should write large enough matrices and pass built-in-panel options to R."""

    calls = []

    def fake_run_r_tool(tool: str, payload: dict[str, object], **kwargs):
        """Capture the xCell payload and return fake score output."""

        calls.append((tool, payload, kwargs))
        written = pd.read_csv(str(payload["mixture_path"]), sep="\t", index_col=0)
        assert written.shape == (5000, 1)
        result = pd.DataFrame({"Epithelial cells": [0.8]}, index=["bulk_1"])
        return result, payload

    monkeypatch.setattr("foli_decon.tools.bulk_tools.run_r_tool", fake_run_r_tool)

    mixture = pd.DataFrame({"bulk_1": [1.0] * 5000}, index=[f"Gene{i}" for i in range(5000)])
    result = run_deconvolution(tool="xcell", mixture=mixture, transform="raw", min_genes=5000)

    assert calls[0][0] == "xcell"
    assert calls[0][1]["arrays"] is False
    assert result.score.loc["bulk_1", "Epithelial cells"] == 0.8
    assert result.metadata["min_genes"] == 5000


def test_run_quantiseq_payload(monkeypatch) -> None:
    """quanTIseq should pass tumor, array, and mRNA scaling options to R."""

    calls = []

    def fake_run_r_tool(tool: str, payload: dict[str, object], **kwargs):
        """Capture the R payload and return a fake fraction table."""

        calls.append((tool, payload, kwargs))
        written = pd.read_csv(str(payload["mixture_path"]), sep="\t", index_col=0)
        assert written.index.tolist() == ["Gene-1", "Gene2"]
        result = pd.DataFrame({"B cells": [0.6], "T cells": [0.4]}, index=["bulk_1"])
        return result, payload

    monkeypatch.setattr("foli_decon.tools.bulk_tools.run_r_tool", fake_run_r_tool)

    mixture = pd.DataFrame({"bulk_1": [10.0, 30.0]}, index=["Gene.1", "Gene2"])
    result = run_deconvolution(
        tool="quantiseq",
        mixture=mixture,
        transform="raw",
        tumor=True,
        arrays=False,
        scale_mrna=False,
    )

    assert calls[0][0] == "quantiseq"
    assert calls[0][1]["tumor"] is True
    assert calls[0][1]["scale_mrna"] is False
    assert result.proportion.loc["bulk_1", "B cells"] == 0.6
    assert result.metadata["transform"] == "raw"


def test_run_mcp_counter_alias_and_payload(monkeypatch) -> None:
    """MCP-counter should accept hyphenated dispatcher names and feature types."""

    calls = []

    def fake_run_r_tool(tool: str, payload: dict[str, object], **kwargs):
        """Capture the R payload and return a fake score table."""

        calls.append((tool, payload, kwargs))
        result = pd.DataFrame({"T cells": [12.0], "B lineage": [8.0]}, index=["bulk_1"])
        return result, payload

    monkeypatch.setattr("foli_decon.tools.bulk_tools.run_r_tool", fake_run_r_tool)

    mixture = pd.DataFrame({"bulk_1": [5.0, 15.0]}, index=["Gene1", "Gene2"])
    result = run_deconvolution(
        tool="MCP-counter",
        mixture=mixture,
        transform="raw",
        feature_types="ENTREZ_ID",
    )

    assert calls[0][0] == "mcp_counter"
    assert calls[0][1]["feature_types"] == "ENTREZ_ID"
    assert calls[0][1]["probesets_path"] == ""
    assert calls[0][1]["genes_path"] == ""
    assert result.tool == "mcp_counter"
    assert result.score.loc["bulk_1", "T cells"] == 12.0


def test_run_psea_marker_signal_payload(monkeypatch) -> None:
    """PSEA should write marker groups and return score-like reference signals."""

    calls = []

    def fake_run_r_tool(tool: str, payload: dict[str, object], **kwargs):
        """Capture the PSEA payload and return fake marker reference signals."""

        calls.append((tool, payload, kwargs))
        markers = pd.read_csv(str(payload["marker_sets_path"]), sep="\t")
        assert markers["cell_type"].tolist() == ["A", "B", "B"]
        assert markers["gene"].tolist() == ["Gene-1", "Gene2", "Gene3"]
        result = pd.DataFrame({"A": [1.2], "B": [0.8]}, index=["bulk_1"])
        return result, payload

    monkeypatch.setattr("foli_decon.tools.bulk_tools.run_r_tool", fake_run_r_tool)

    mixture = pd.DataFrame({"bulk_1": [10.0, 30.0, 20.0]}, index=["Gene.1", "Gene2", "Gene3"])
    result = run_deconvolution(
        tool="psea",
        mixture=mixture,
        marker_sets={"A": ["Gene.1"], "B": [["Gene2", "Gene3"]]},
        transform="raw",
        sample_subset=["bulk_1"],
        target_mean=100.0,
    )

    assert calls[0][0] == "psea"
    assert calls[0][1]["sample_subset"] == ["bulk_1"]
    assert calls[0][1]["target_mean"] == 100.0
    assert result.score.loc["bulk_1", "A"] == 1.2
    assert result.metadata["output_semantics"] == "PSEA marker reference signal; not a compositional fraction"


def test_run_cellcode_payload_and_score(monkeypatch, tmp_path) -> None:
    """CellCODE should return surrogate proportion variables as score output."""

    calls = []

    def fake_run_r_tool(tool: str, payload: dict[str, object], **kwargs):
        """Capture CellCODE payload and return fake SPV score output."""

        calls.append((tool, payload, kwargs))
        mixture = pd.read_csv(str(payload["mixture_path"]), sep="\t", index_col=0)
        signature = pd.read_csv(str(payload["signature_path"]), sep="\t", index_col=0)
        groups = pd.read_csv(str(payload["groups_path"]), sep="\t")
        assert mixture.index.tolist() == ["Gene-1", "Gene2"]
        assert signature.index.tolist() == ["Gene-1", "Gene2"]
        assert groups["group"].tolist() == ["control", "case"]
        assert payload["method"] == "mixed"
        assert payload["max_markers"] == 2
        result = pd.DataFrame({"A": [1.1, -0.4], "B": [0.2, 0.8]}, index=["bulk_1", "bulk_2"])
        return result, payload

    monkeypatch.setattr("foli_decon.tools.bulk_tools.run_r_tool", fake_run_r_tool)

    mixture = pd.DataFrame(
        {
            "bulk_1": [10.0, 30.0, 50.0],
            "bulk_2": [20.0, 40.0, 60.0],
        },
        index=["Gene.1", "Gene2", "Gene3"],
    )
    signature = pd.DataFrame(
        {
            "A": [1.0, 2.0],
            "B": [3.0, 4.0],
        },
        index=["Gene.1", "Gene2"],
    )
    groups = pd.Series(["control", "case"], index=mixture.columns)

    result = run_deconvolution(
        tool="cellcode",
        mixture=mixture,
        signature=signature,
        groups=groups,
        max_markers=2,
        r_env="foli-decon-r44",
        r_output_path=str(tmp_path / "cellcode.tsv"),
    )

    assert calls[0][0] == "cellcode"
    assert calls[0][2]["rscript_command"] == ["mamba", "run", "-n", "foli-decon-r44", "Rscript"]
    assert result.proportion is None
    assert result.score.loc["bulk_2", "B"] == 0.8
    assert result.output_paths["score"] == str(tmp_path / "cellcode.tsv")


def test_run_abis_payload(monkeypatch, tmp_path) -> None:
    """ABIS should stage the native signature matrix and return fraction-scaled output."""

    calls = []

    def fake_run_r_tool(tool: str, payload: dict[str, object], **kwargs):
        """Capture the ABIS payload and return fake native percent output."""

        calls.append((tool, payload, kwargs))
        mixture = pd.read_csv(str(payload["mixture_path"]), sep="\t", index_col=0)
        signature = pd.read_csv(str(payload["signature_path"]), sep="\t", index_col=0)
        assert mixture.index.tolist() == ["Gene-1", "Gene2"]
        assert signature.index.tolist() == ["Gene-1", "Gene2"]
        assert payload["technology"] == "rnaseq"
        assert payload["target_path"] == ""
        assert payload["maxit"] == 100
        result = pd.DataFrame({"Monocytes C": [10.0], "NK": [20.0]}, index=["bulk_1"])
        return result, payload

    monkeypatch.setattr("foli_decon.tools.bulk_tools.run_r_tool", fake_run_r_tool)

    mixture = pd.DataFrame({"bulk_1": [10.0, 30.0]}, index=["Gene.1", "Gene2"])
    signature = pd.DataFrame(
        {
            "Monocytes C": [1.0, 2.0],
            "NK": [3.0, 4.0],
        },
        index=["Gene.1", "Gene2"],
    )

    result = run_deconvolution(
        tool="abis",
        mixture=mixture,
        signature=signature,
        transform="raw",
        technology="rnaseq",
        r_env="foli-decon",
        output_path=str(tmp_path / "abis.tsv"),
        log_path=str(tmp_path / "abis.log"),
    )

    assert calls[0][0] == "abis"
    assert calls[0][2]["rscript_command"] == ["mamba", "run", "-n", "foli-decon", "Rscript"]
    assert calls[0][2]["output_path"] == str(tmp_path / "abis.tsv")
    assert result.proportion.loc["bulk_1", "Monocytes C"] == 0.1
    assert result.metadata["native_output_scale"] == "percent"
    assert result.metadata["output_scale"] == "fraction"
    assert result.output_paths["proportion"] == str(tmp_path / "abis.tsv")


def test_run_estimate_payload(monkeypatch, tmp_path) -> None:
    """ESTIMATE should return scores instead of proportions."""

    calls = []

    def fake_run_r_tool(tool: str, payload: dict[str, object], **kwargs):
        """Capture the ESTIMATE payload and return fake score output."""

        calls.append((tool, payload, kwargs))
        mixture = pd.read_csv(str(payload["mixture_path"]), sep="\t", index_col=0)
        assert mixture.index.tolist() == ["Gene-1", "Gene2"]
        assert payload["platform"] == "illumina"
        result = pd.DataFrame(
            {"StromalScore": [1.0], "ImmuneScore": [2.0], "ESTIMATEScore": [3.0]},
            index=["bulk_1"],
        )
        return result, payload

    monkeypatch.setattr("foli_decon.tools.bulk_tools.run_r_tool", fake_run_r_tool)

    mixture = pd.DataFrame({"bulk_1": [10.0, 30.0]}, index=["Gene.1", "Gene2"])
    result = run_deconvolution(
        tool="estimate",
        mixture=mixture,
        transform="raw",
        platform="illumina",
        output_path=str(tmp_path / "estimate.tsv"),
    )

    assert calls[0][0] == "estimate"
    assert calls[0][2]["output_path"] == str(tmp_path / "estimate.tsv")
    assert result.proportion is None
    assert result.score.loc["bulk_1", "ESTIMATEScore"] == 3.0
    assert result.metadata["output_semantics"].startswith("ESTIMATE stromal")
    assert result.output_paths["score"] == str(tmp_path / "estimate.tsv")


def test_run_autogenes_sidecar_payload(monkeypatch, tmp_path) -> None:
    """AutoGeneS should align genes and call the Python sidecar runner."""

    calls = []

    def fake_run_command(command: list[str], **kwargs) -> None:
        """Capture the AutoGeneS command and write fake sidecar outputs."""

        calls.append((command, kwargs))
        payload = json.loads(Path(command[-1]).read_text(encoding="utf-8"))
        mixture = pd.read_csv(payload["mixture_path"], sep="\t", index_col=0)
        signature = pd.read_csv(payload["signature_path"], sep="\t", index_col=0)
        assert mixture.index.tolist() == ["Gene-1", "Gene2"]
        assert signature.index.tolist() == ["Gene-1", "Gene2"]
        assert payload["ngen"] == 3
        assert payload["nfeatures"] == 2
        pd.DataFrame({"A": [0.75], "B": [0.25]}, index=["bulk_1"]).to_csv(
            payload["output_path"],
            sep="\t",
        )
        pd.DataFrame({"gene": ["Gene-1", "Gene2"]}).to_csv(
            payload["selected_features_path"],
            sep="\t",
            index=False,
        )

    monkeypatch.setattr("foli_decon.tools.bulk_tools.run_command", fake_run_command)

    mixture = pd.DataFrame({"bulk_1": [10.0, 20.0, 30.0]}, index=["Gene.1", "Gene2", "Gene3"])
    signature = pd.DataFrame(
        {
            "A": [1.0, 2.0, 3.0],
            "B": [3.0, 2.0, 1.0],
        },
        index=["Gene.1", "Gene2", "Gene4"],
    )

    result = run_deconvolution(
        tool="autogenes",
        mixture=mixture,
        signature=signature,
        mixture_transform="raw",
        signature_transform="raw",
        python_env="foli-decon-py311-cpu",
        ngen=3,
        nfeatures=2,
        selected_features_path=str(tmp_path / "selected.tsv"),
    )

    assert calls[0][0][:5] == ["mamba", "run", "-n", "foli-decon-py311-cpu", "python"]
    assert result.proportion.loc["bulk_1", "A"] == 0.75
    assert result.selected_features["gene"].tolist() == ["Gene-1", "Gene2"]
    assert result.output_paths["selected_features"] == str(tmp_path / "selected.tsv")


def test_select_features_autogenes_sidecar_payload(monkeypatch, tmp_path) -> None:
    """AutoGeneS feature selection should be available without deconvolution."""

    calls = []

    def fake_run_command(command: list[str], **kwargs) -> None:
        """Capture the AutoGeneS selection command and write fake selected genes."""

        calls.append((command, kwargs))
        payload = json.loads(Path(command[-1]).read_text(encoding="utf-8"))
        signature = pd.read_csv(payload["signature_path"], sep="\t", index_col=0)
        assert signature.index.tolist() == ["Gene-1", "Gene2", "Gene4"]
        assert payload["nfeatures"] == 2
        pd.DataFrame({"gene": ["Gene-1", "Gene2"]}).to_csv(
            payload["selected_features_path"],
            sep="\t",
            index=False,
        )

    monkeypatch.setattr("foli_decon.tools.bulk_tools.run_command", fake_run_command)

    signature = pd.DataFrame(
        {
            "A": [1.0, 2.0, 3.0],
            "B": [3.0, 2.0, 1.0],
        },
        index=["Gene.1", "Gene2", "Gene4"],
    )

    direct_result = select_features_autogenes(
        signature=signature,
        signature_transform="raw",
        python_env="foli-decon-py311-cpu",
        nfeatures=2,
        selected_features_path=str(tmp_path / "selected.tsv"),
    )
    dispatched_result = select_features(
        tool="autogenes",
        signature=signature,
        signature_transform="raw",
        nfeatures=2,
    )

    assert supported_feature_selectors() == ("autogenes", "markermap", "scgenefit")
    assert calls[0][0][:5] == ["mamba", "run", "-n", "foli-decon-py311-cpu", "python"]
    assert direct_result.selected_genes.tolist() == ["Gene-1", "Gene2"]
    assert direct_result.output_paths["selected_features"] == str(tmp_path / "selected.tsv")
    assert dispatched_result.selected_genes.tolist() == ["Gene-1", "Gene2"]


def test_select_features_scgenefit_sidecar_payload(monkeypatch, tmp_path) -> None:
    """scGeneFit feature selection should pass dense cell-label data to its sidecar."""

    calls = []

    def fake_run_command(command: list[str], **kwargs) -> None:
        """Capture the scGeneFit command and write fake selected genes."""

        calls.append((command, kwargs))
        payload = json.loads(Path(command[-1]).read_text(encoding="utf-8"))
        scrna = pd.read_csv(payload["scrna_path"], sep="\t", index_col=0)
        labels = pd.read_csv(payload["labels_path"], sep="\t")
        assert scrna.index.tolist() == ["Gene-1", "Gene2", "Gene3"]
        assert labels["cell_type"].tolist() == ["A", "A", "B", "B"]
        assert payload["nfeatures"] == 2
        assert payload["method"] in {"pairwise_centers", "centers"}
        pd.DataFrame({"gene": ["Gene-1", "Gene3"]}).to_csv(
            payload["selected_features_path"],
            sep="\t",
            index=False,
        )
        pd.DataFrame(
            {
                "rank": [1, 2],
                "gene": ["Gene-1", "Gene3"],
                "feature_index": [0, 2],
            }
        ).to_csv(payload["ranking_path"], sep="\t", index=False)

    monkeypatch.setattr("foli_decon.tools.bulk_tools.run_command", fake_run_command)

    scrna_counts = pd.DataFrame(
        {
            "cell_1": [10.0, 2.0, 1.0],
            "cell_2": [9.0, 1.0, 1.0],
            "cell_3": [1.0, 2.0, 8.0],
            "cell_4": [1.0, 1.0, 9.0],
        },
        index=["Gene.1", "Gene2", "Gene3"],
    )
    cell_types = pd.Series(["A", "A", "B", "B"], index=scrna_counts.columns)

    direct_result = select_features_scgenefit(
        scrna_counts=scrna_counts,
        cell_types=cell_types,
        scrna_transform="raw",
        python_env="foli-decon-py311-cpu",
        nfeatures=2,
        method="pairwise_centers",
        selected_features_path=str(tmp_path / "scgenefit_selected.tsv"),
        ranking_path=str(tmp_path / "scgenefit_ranking.tsv"),
    )
    dispatched_result = select_features(
        tool="scGeneFit",
        scrna_counts=scrna_counts,
        cell_types=cell_types,
        scrna_transform="raw",
        nfeatures=2,
    )

    assert calls[0][0][:5] == ["mamba", "run", "-n", "foli-decon-py311-cpu", "python"]
    assert direct_result.selected_genes.tolist() == ["Gene-1", "Gene3"]
    assert direct_result.ranking["feature_index"].tolist() == [0, 2]
    assert direct_result.output_paths["selected_features"] == str(tmp_path / "scgenefit_selected.tsv")
    assert direct_result.output_paths["ranking"] == str(tmp_path / "scgenefit_ranking.tsv")
    assert dispatched_result.tool == "scgenefit"


def test_select_features_markermap_sidecar_payload(monkeypatch, tmp_path) -> None:
    """MarkerMap feature selection should pass training controls to its sidecar."""

    calls = []

    def fake_run_command(command: list[str], **kwargs) -> None:
        """Capture the MarkerMap command and write fake selected genes."""

        calls.append((command, kwargs))
        payload = json.loads(Path(command[-1]).read_text(encoding="utf-8"))
        scrna = pd.read_csv(payload["scrna_path"], sep="\t", index_col=0)
        labels = pd.read_csv(payload["labels_path"], sep="\t")
        assert scrna.index.tolist() == ["Gene-1", "Gene2", "Gene3", "Gene4"]
        assert labels["cell_type"].tolist() == ["A", "A", "B", "B"]
        assert payload["nfeatures"] == 2
        assert payload["hidden_layer_size"] in {4, 7}
        assert payload["max_epochs"] in {3, 10}
        pd.DataFrame({"gene": ["Gene-1", "Gene4"]}).to_csv(
            payload["selected_features_path"],
            sep="\t",
            index=False,
        )
        pd.DataFrame(
            {
                "rank": [1, 2],
                "gene": ["Gene-1", "Gene4"],
                "feature_index": [0, 3],
            }
        ).to_csv(payload["ranking_path"], sep="\t", index=False)

    monkeypatch.setattr("foli_decon.tools.bulk_tools.run_command", fake_run_command)

    scrna_counts = pd.DataFrame(
        {
            "cell_1": [10.0, 2.0, 1.0, 3.0],
            "cell_2": [9.0, 1.0, 1.0, 2.0],
            "cell_3": [1.0, 2.0, 8.0, 6.0],
            "cell_4": [1.0, 1.0, 9.0, 7.0],
        },
        index=["Gene.1", "Gene2", "Gene3", "Gene4"],
    )
    cell_types = pd.Series(["A", "A", "B", "B"], index=scrna_counts.columns)

    direct_result = select_features_markermap(
        scrna_counts=scrna_counts,
        cell_types=cell_types,
        scrna_transform="raw",
        python_env="foli-decon-py311-cpu",
        nfeatures=2,
        hidden_layer_size=7,
        min_epochs=1,
        max_epochs=3,
        selected_features_path=str(tmp_path / "markermap_selected.tsv"),
        ranking_path=str(tmp_path / "markermap_ranking.tsv"),
    )
    dispatched_result = select_features(
        tool="markermap",
        scrna_counts=scrna_counts,
        cell_types=cell_types,
        scrna_transform="raw",
        nfeatures=2,
    )

    assert calls[0][0][:5] == ["mamba", "run", "-n", "foli-decon-py311-cpu", "python"]
    assert direct_result.selected_genes.tolist() == ["Gene-1", "Gene4"]
    assert direct_result.ranking["feature_index"].tolist() == [0, 3]
    assert direct_result.output_paths["selected_features"] == str(tmp_path / "markermap_selected.tsv")
    assert direct_result.output_paths["ranking"] == str(tmp_path / "markermap_ranking.tsv")
    assert dispatched_result.tool == "markermap"


def test_run_blade_signature_sidecar_payload(monkeypatch, tmp_path) -> None:
    """BLADE should write aligned mean/sd references and call the Python sidecar."""

    calls = []

    def fake_run_command(command: list[str], **kwargs) -> None:
        """Capture the BLADE payload and write fake sidecar outputs."""

        calls.append((command, kwargs))
        payload = json.loads(Path(command[-1]).read_text(encoding="utf-8"))
        mixture = pd.read_csv(payload["mixture_path"], sep="\t", index_col=0)
        signature_mean = pd.read_csv(payload["signature_mean_path"], sep="\t", index_col=0)
        signature_sd = pd.read_csv(payload["signature_sd_path"], sep="\t", index_col=0)
        assert mixture.index.tolist() == ["Gene-1", "Gene2"]
        assert signature_mean.index.tolist() == ["Gene-1", "Gene2"]
        assert signature_sd.index.tolist() == ["Gene-1", "Gene2"]
        assert payload["alphas"] == [1.0]
        assert payload["n_rep"] == 1
        pd.DataFrame({"A": [0.8], "B": [0.2]}, index=["bulk_1"]).to_csv(
            payload["output_path"],
            sep="\t",
        )

    monkeypatch.setattr("foli_decon.tools.bulk_tools.run_command", fake_run_command)

    mixture = pd.DataFrame({"bulk_1": [10.0, 20.0, 30.0]}, index=["Gene.1", "Gene2", "Gene3"])
    signature = pd.DataFrame({"A": [1.0, 2.0], "B": [3.0, 4.0]}, index=["Gene.1", "Gene2"])
    signature_sd = pd.DataFrame({"A": [0.1, 0.2], "B": [0.3, 0.4]}, index=["Gene.1", "Gene2"])

    result = run_deconvolution(
        tool="blade",
        mixture=mixture,
        signature=signature,
        signature_sd=signature_sd,
        mixture_transform="raw",
        signature_transform="raw",
        python_env="foli-decon-py311-cpu",
        alphas=[1.0],
        alpha0s=[0.1],
        kappa0s=[1.0],
        sigma_ys=[0.5],
        n_rep=1,
        n_rep_final=1,
        n_jobs=1,
        extra_output_dir=str(tmp_path / "blade_extra"),
    )

    assert calls[0][0][:5] == ["mamba", "run", "-n", "foli-decon-py311-cpu", "python"]
    assert result.proportion.loc["bulk_1", "A"] == 0.8
    assert result.output_paths["extra_output_dir"] == str(tmp_path / "blade_extra")


def test_run_consensus_tme_natural_alias(monkeypatch) -> None:
    """ConsensusTME should dispatch when passed without an underscore."""

    calls = []

    def fake_run_r_tool(tool: str, payload: dict[str, object], **kwargs):
        """Capture the ConsensusTME payload and return fake score output."""

        calls.append((tool, payload, kwargs))
        result = pd.DataFrame({"T cells": [3.0]}, index=["bulk_1"])
        return result, payload

    monkeypatch.setattr("foli_decon.tools.bulk_tools.run_r_tool", fake_run_r_tool)

    mixture = pd.DataFrame({"bulk_1": [5.0, 15.0]}, index=["Gene1", "Gene2"])
    result = run_deconvolution(tool="ConsensusTME", mixture=mixture, transform="raw")

    assert calls[0][0] == "consensus_tme"
    assert result.score.loc["bulk_1", "T cells"] == 3.0


def test_run_cdseq_payload(monkeypatch) -> None:
    """CDSeq should build a pseudobulk reference and pass sidecar R env settings."""

    calls = []

    def fake_run_r_tool(tool: str, payload: dict[str, object], **kwargs):
        """Capture CDSeq payload and return fake latent/reference-labeled proportions."""

        calls.append((tool, payload, kwargs))
        mixture = pd.read_csv(str(payload["mixture_path"]), sep="\t", index_col=0)
        reference = pd.read_csv(str(payload["reference_gep_path"]), sep="\t", index_col=0)
        assert mixture.index.tolist() == ["Gene-1", "Gene2"]
        assert reference.columns.tolist() == ["A", "B"]
        result = pd.DataFrame({"A": [0.65], "B": [0.35]}, index=["bulk_1"])
        return result, payload

    monkeypatch.setattr("foli_decon.tools.bulk_tools.run_r_tool", fake_run_r_tool)

    mixture = pd.DataFrame({"bulk_1": [10.0, 30.0]}, index=["Gene.1", "Gene2"])
    scrna = pd.DataFrame(
        {
            "cell_1": [1.0, 2.0],
            "cell_2": [3.0, 4.0],
        },
        index=["Gene.1", "Gene2"],
    )
    cell_types = pd.Series(["A", "B"], index=scrna.columns)

    result = run_deconvolution(
        tool="cdseq",
        mixture=mixture,
        scrna_counts=scrna,
        cell_types=cell_types,
        mcmc_iterations=5,
        r_env="foli-decon-r44",
    )

    assert calls[0][0] == "cdseq"
    assert calls[0][1]["cell_type_number"] == [2]
    assert calls[0][1]["mcmc_iterations"] == 5
    assert calls[0][2]["rscript_command"] == ["mamba", "run", "-n", "foli-decon-r44", "Rscript"]
    assert result.proportion.loc["bulk_1", "A"] == 0.65


def test_run_music2_matrix_payload(monkeypatch) -> None:
    """MuSiC2 matrix mode should align control, case, and scRNA genes before R."""

    calls = []

    def fake_run_r_tool(tool: str, payload: dict[str, object], **kwargs):
        """Capture matrix-mode payload paths while temp files still exist."""

        calls.append((tool, payload, kwargs))
        control = pd.read_csv(str(payload["control_mixture_path"]), sep="\t", index_col=0)
        case = pd.read_csv(str(payload["case_mixture_path"]), sep="\t", index_col=0)
        scrna = pd.read_csv(str(payload["scrna_path"]), sep="\t", index_col=0)
        assert control.index.tolist() == ["Gene-1", "Gene2"]
        assert case.index.tolist() == ["Gene-1", "Gene2"]
        assert scrna.index.tolist() == ["Gene-1", "Gene2"]
        result = pd.DataFrame({"A": [0.7], "B": [0.3]}, index=["case_1"])
        return result, payload

    monkeypatch.setattr("foli_decon.tools.scrna_tools.run_r_tool", fake_run_r_tool)

    control = pd.DataFrame({"ctrl_1": [10.0, 20.0, 30.0]}, index=["Gene.1", "Gene2", "Gene3"])
    case = pd.DataFrame({"case_1": [15.0, 25.0]}, index=["Gene.1", "Gene2"])
    scrna = pd.DataFrame(
        {
            "cell_1": [1.0, 2.0, 3.0],
            "cell_2": [2.0, 1.0, 4.0],
        },
        index=["Gene.1", "Gene2", "Gene4"],
    )
    cell_types = pd.Series(["A", "B"], index=scrna.columns)
    batch_ids = pd.Series(["d1", "d2"], index=scrna.columns)

    result = run_music2(
        control_mixture=control,
        case_mixture=case,
        scrna_counts=scrna,
        cell_types=cell_types,
        batch_ids=batch_ids,
        n_resample=3,
        cutoff_fc=1.5,
    )

    assert calls[0][0] == "music2"
    assert calls[0][1]["n_resample"] == 3
    assert calls[0][1]["cutoff_fc"] == 1.5
    assert result.metadata["reference_mode"] == "matrix"
    assert result.proportion.loc["case_1", "A"] == 0.7


def test_run_mead_payload_and_uncertainty(monkeypatch, tmp_path) -> None:
    """MEAD should build SingleCellExperiment inputs and read p_hat standard errors."""

    calls = []

    def fake_run_r_tool(tool: str, payload: dict[str, object], **kwargs):
        """Capture MEAD payload paths and write fake uncertainty output."""

        calls.append((tool, payload, kwargs))
        mixture = pd.read_csv(str(payload["mixture_path"]), sep="\t", index_col=0)
        scrna = pd.read_csv(str(payload["scrna_path"]), sep="\t", index_col=0)
        groups = pd.read_csv(str(payload["groups_path"]), sep="\t")
        assert mixture.index.tolist() == ["Gene-1", "Gene2"]
        assert scrna.index.tolist() == ["Gene-1", "Gene2"]
        assert groups["group"].tolist() == ["case"]
        assert payload["filter_gene"] is False
        assert payload["gene_thresh"] == 0.0
        assert payload["calc_var"] is True
        pd.DataFrame({"A": [0.01], "B": [0.02]}, index=["bulk_1"]).to_csv(
            str(payload["uncertainty_path"]),
            sep="\t",
        )
        result = pd.DataFrame({"A": [0.7], "B": [0.3]}, index=["bulk_1"])
        return result, payload

    monkeypatch.setattr("foli_decon.tools.scrna_tools.run_r_tool", fake_run_r_tool)

    mixture = pd.DataFrame({"bulk_1": [10.0, 20.0, 30.0]}, index=["Gene.1", "Gene2", "Gene3"])
    scrna = pd.DataFrame(
        {
            "cell_1": [1.0, 2.0, 3.0],
            "cell_2": [2.0, 1.0, 4.0],
        },
        index=["Gene.1", "Gene2", "Gene4"],
    )
    cell_types = pd.Series(["A", "B"], index=scrna.columns)
    batch_ids = pd.Series(["d1", "d2"], index=scrna.columns)
    groups = pd.Series(["case"], index=mixture.columns)

    result = run_deconvolution(
        tool="mead",
        mixture=mixture,
        scrna_counts=scrna,
        cell_types=cell_types,
        batch_ids=batch_ids,
        groups=groups,
        r_env="foli-decon-r44",
        r_output_path=str(tmp_path / "mead.tsv"),
        extra_output_dir=str(tmp_path / "mead_extra"),
    )

    assert calls[0][0] == "mead"
    assert calls[0][2]["rscript_command"] == ["mamba", "run", "-n", "foli-decon-r44", "Rscript"]
    assert result.proportion.loc["bulk_1", "A"] == 0.7
    assert result.uncertainty.loc["bulk_1", "B"] == 0.02
    assert result.output_paths["proportion"] == str(tmp_path / "mead.tsv")
    assert result.output_paths["extra_output_dir"] == str(tmp_path / "mead_extra")


def test_run_music2_sce_payload(monkeypatch, tmp_path) -> None:
    """MuSiC2 SCE mode should pass metadata columns and sidecar R env through."""

    calls = []

    def fake_run_r_tool(tool: str, payload: dict[str, object], **kwargs):
        """Capture SCE-mode payload and return a fake fraction table."""

        calls.append((tool, payload, kwargs))
        result = pd.DataFrame({"A": [0.2], "B": [0.8]}, index=["case_1"])
        return result, payload

    monkeypatch.setattr("foli_decon.tools.scrna_tools.run_r_tool", fake_run_r_tool)

    control = pd.DataFrame({"ctrl_1": [10.0, 20.0]}, index=["Gene1", "Gene2"])
    case = pd.DataFrame({"case_1": [15.0, 25.0]}, index=["Gene1", "Gene2"])
    sce_path = tmp_path / "reference.rds"

    result = run_music2(
        control_mixture=control,
        case_mixture=case,
        scrna_sce_rds_path=str(sce_path),
        cell_type_column="celltype",
        batch_id_column="donor",
        r_env="foli-decon-r44",
        extra_output_dir=str(tmp_path / "extras"),
    )

    assert calls[0][0] == "music2_sce"
    assert calls[0][1]["scrna_sce_rds_path"] == str(sce_path)
    assert calls[0][1]["cell_type_column"] == "celltype"
    assert calls[0][1]["batch_id_column"] == "donor"
    assert calls[0][2]["rscript_command"] == ["mamba", "run", "-n", "foli-decon-r44", "Rscript"]
    assert result.metadata["reference_mode"] == "sce_rds"
    assert result.proportion.loc["case_1", "B"] == 0.8


def test_run_scdc_payload(monkeypatch) -> None:
    """SCDC should pass donor-aware ExpressionSet metadata through R payload files."""

    calls = []

    def fake_run_r_tool(tool: str, payload: dict[str, object], **kwargs):
        """Capture SCDC payload paths and return fake proportions."""

        calls.append((tool, payload, kwargs))
        scrna = pd.read_csv(str(payload["scrna_path"]), sep="\t", index_col=0)
        assert scrna.index.tolist() == ["Gene1", "Gene2"]
        result = pd.DataFrame({"A": [0.55], "B": [0.45]}, index=["bulk_1"])
        return result, payload

    monkeypatch.setattr("foli_decon.tools.scrna_tools.run_r_tool", fake_run_r_tool)

    mixture = pd.DataFrame({"bulk_1": [10.0, 20.0]}, index=["Gene1", "Gene2"])
    scrna = pd.DataFrame(
        {
            "cell_1": [1.0, 2.0],
            "cell_2": [2.0, 1.0],
        },
        index=["Gene1", "Gene2"],
    )
    cell_types = pd.Series(["A", "B"], index=scrna.columns)
    batch_ids = pd.Series(["d1", "d2"], index=scrna.columns)

    result = run_deconvolution(
        tool="scdc",
        mixture=mixture,
        scrna_counts=scrna,
        cell_types=cell_types,
        batch_ids=batch_ids,
        ct_sub=["A", "B"],
        iter_max=5,
    )

    assert calls[0][0] == "scdc"
    assert calls[0][1]["ct_sub"] == ["A", "B"]
    assert calls[0][1]["iter_max"] == 5
    assert result.proportion.loc["bulk_1", "A"] == 0.55


def test_run_instaprism_alias_payload(monkeypatch) -> None:
    """InstaPrism should accept the instaprime alias and optional cell states."""

    calls = []

    def fake_run_r_tool(tool: str, payload: dict[str, object], **kwargs):
        """Capture InstaPrism payload and return fake proportions."""

        calls.append((tool, payload, kwargs))
        cell_states = pd.read_csv(str(payload["cell_states_path"]), sep="\t")
        assert cell_states.iloc[:, 0].tolist() == ["A_1", "B_1"]
        result = pd.DataFrame({"A": [0.25], "B": [0.75]}, index=["bulk_1"])
        return result, payload

    monkeypatch.setattr("foli_decon.tools.scrna_tools.run_r_tool", fake_run_r_tool)

    mixture = pd.DataFrame({"bulk_1": [10.0, 20.0]}, index=["Gene1", "Gene2"])
    scrna = pd.DataFrame(
        {
            "cell_1": [1.0, 2.0],
            "cell_2": [2.0, 1.0],
        },
        index=["Gene1", "Gene2"],
    )
    cell_types = pd.Series(["A", "B"], index=scrna.columns)
    cell_states = pd.Series(["A_1", "B_1"], index=scrna.columns)

    result = run_deconvolution(
        tool="instaprime",
        mixture=mixture,
        scrna_counts=scrna,
        cell_types=cell_types,
        cell_states=cell_states,
        r_env="foli-decon-r44",
        write_z=True,
    )

    assert calls[0][0] == "instaprism"
    assert calls[0][1]["write_z"] is True
    assert calls[0][2]["rscript_command"] == ["mamba", "run", "-n", "foli-decon-r44", "Rscript"]
    assert result.metadata["cell_state_mode"] == "provided"
    assert result.proportion.loc["bulk_1", "B"] == 0.75


def test_run_tape_sidecar_payload(monkeypatch) -> None:
    """TAPE should write cell-type-indexed scRNA references and sidecar payloads."""

    calls = []

    def fake_run_command(command: list[str], **kwargs) -> None:
        """Capture the TAPE payload and write fake fraction output."""

        calls.append((command, kwargs))
        payload = json.loads(Path(command[-1]).read_text(encoding="utf-8"))
        mixture = pd.read_csv(payload["mixture_path"], sep="\t", index_col=0)
        scrna = pd.read_csv(payload["scrna_path"], sep="\t", index_col=0)
        assert mixture.index.tolist() == ["bulk_1"]
        assert mixture.columns.tolist() == ["Gene1", "Gene2"]
        assert scrna.index.tolist() == ["A", "B"]
        assert scrna.columns.tolist() == ["Gene1", "Gene2"]
        assert payload["epochs"] == 2
        pd.DataFrame({"A": [0.35], "B": [0.65]}, index=["bulk_1"]).to_csv(
            payload["output_path"],
            sep="\t",
        )

    monkeypatch.setattr("foli_decon.tools.scrna_tools.run_command", fake_run_command)

    mixture = pd.DataFrame({"bulk_1": [10.0, 20.0]}, index=["Gene1", "Gene2"])
    scrna = pd.DataFrame(
        {
            "cell_1": [1.0, 2.0],
            "cell_2": [2.0, 1.0],
        },
        index=["Gene1", "Gene2"],
    )
    cell_types = pd.Series(["A", "B"], index=scrna.columns)

    result = run_deconvolution(
        tool="tape",
        mixture=mixture,
        scrna_counts=scrna,
        cell_types=cell_types,
        python_env="foli-decon-py311-cpu",
        epochs=2,
        batch_size=4,
    )

    assert calls[0][0][:5] == ["mamba", "run", "-n", "foli-decon-py311-cpu", "python"]
    assert result.proportion.loc["bulk_1", "B"] == 0.65
    assert result.metadata["python_env"] == "foli-decon-py311-cpu"


def test_run_scaden_sidecar_commands(monkeypatch) -> None:
    """Scaden should write native text inputs and call simulate/process/train/predict."""

    calls = []

    def fake_run_command(command: list[str], cwd: str | None = None, log_path: str | None = None):
        """Capture Scaden commands and create fake prediction output."""

        calls.append(command)
        if "simulate" in command:
            data_dir = Path(command[command.index("--data") + 1])
            counts = pd.read_csv(data_dir / "reference_counts.txt", sep="\t", index_col=0)
            celltypes = pd.read_csv(data_dir / "reference_celltypes.txt", sep="\t")
            assert counts.index.tolist() == [0, 1]
            assert counts.columns.tolist() == ["Gene1", "Gene2"]
            assert celltypes["Celltype"].tolist() == ["A", "B"]
        if "predict" in command:
            outname = Path(command[command.index("--outname") + 1])
            pd.DataFrame({"A": [0.4], "B": [0.6]}, index=["bulk_1"]).to_csv(outname, sep="\t")

    monkeypatch.setattr("foli_decon.tools.scrna_tools.run_command", fake_run_command)

    mixture = pd.DataFrame({"bulk_1": [10.0, 20.0]}, index=["Gene1", "Gene2"])
    scrna = pd.DataFrame(
        {
            "cell_1": [1.0, 2.0],
            "cell_2": [2.0, 1.0],
        },
        index=["Gene1", "Gene2"],
    )
    cell_types = pd.Series(["A", "B"], index=scrna.columns)

    result = run_deconvolution(
        tool="scaden",
        mixture=mixture,
        scrna_counts=scrna,
        cell_types=cell_types,
        python_env="foli-decon-deep-cpu",
        n_training_samples=6,
        cells_per_sample=4,
        train_steps=2,
    )

    assert [command[5] for command in calls] == ["simulate", "process", "train", "predict"]
    assert result.proportion.loc["bulk_1", "B"] == 0.6
    assert result.metadata["python_env"] == "foli-decon-deep-cpu"


def test_run_dissect_sidecar_payload(monkeypatch, tmp_path) -> None:
    """DISSECT should write a h5ad payload and parse native fraction output."""

    calls = []

    def fake_run_command(command: list[str], cwd: str | None = None, log_path: str | None = None):
        """Capture DISSECT sidecar payload and create fake native outputs."""

        calls.append(command)
        payload = json.loads(Path(command[-1]).read_text())
        assert payload["n_models"] == 1
        assert payload["train_steps"] == 2
        assert payload["batch_col"] == "batch"
        scrna = pd.read_csv(payload["scrna_matrix_path"], sep="\t", index_col=0)
        metadata = pd.read_csv(payload["scrna_metadata_path"], sep="\t")
        assert scrna.index.tolist() == ["Gene1", "Gene2"]
        assert metadata["cell_id"].tolist() == ["cell_1", "cell_2"]
        assert metadata["batch"].tolist() == ["d1", "d2"]
        experiment_folder = Path(payload["experiment_folder"])
        experiment_folder.mkdir(parents=True)
        pd.DataFrame({"A": [0.7], "B": [0.3]}, index=["bulk_1"]).to_csv(
            experiment_folder / "dissect_fractions.txt",
            sep="\t",
        )
        pd.DataFrame({"A": [1.2], "B": [-0.2]}, index=["bulk_1"]).to_csv(
            experiment_folder / "dissect_scores.txt",
            sep="\t",
        )
        (experiment_folder / "main_config.json").write_text("{}", encoding="utf-8")

    monkeypatch.setattr("foli_decon.tools.scrna_tools.run_command", fake_run_command)

    mixture = pd.DataFrame({"bulk_1": [10.0, 20.0]}, index=["Gene1", "Gene2"])
    scrna = pd.DataFrame(
        {
            "cell_1": [1.0, 2.0],
            "cell_2": [2.0, 1.0],
        },
        index=["Gene1", "Gene2"],
    )
    cell_types = pd.Series(["A", "B"], index=scrna.columns)
    batch_ids = pd.Series(["d1", "d2"], index=scrna.columns)

    result = run_deconvolution(
        tool="dissect",
        mixture=mixture,
        scrna_counts=scrna,
        cell_types=cell_types,
        batch_ids=batch_ids,
        python_env="foli-decon-deep-cpu",
        n_models=1,
        train_steps=2,
        n_training_samples=6,
        cells_per_sample=4,
        extra_output_dir=str(tmp_path / "dissect_extra"),
    )

    assert calls[0][:5] == ["mamba", "run", "-n", "foli-decon-deep-cpu", "python"]
    assert result.proportion.loc["bulk_1", "A"] == 0.7
    assert (tmp_path / "dissect_extra" / "dissect_scores.txt").read_text()


def test_run_blue_sidecar_payload(monkeypatch, tmp_path) -> None:
    """BLUE should stage a pipeline payload and parse staged native outputs."""

    calls = []

    def fake_run_command(command: list[str], cwd: str | None = None, log_path: str | None = None):
        """Capture BLUE sidecar payload and create fake native outputs."""

        calls.append(command)
        payload = json.loads(Path(command[-1]).read_text())
        assert payload["epochs"] == 1
        assert payload["samplenum_all_train"] == 8
        assert payload["celltype_mapping"] == {"A": ["A"], "B": ["B"]}
        assert Path(payload["blue_repo_path"]).name == "BLUE"
        scrna = pd.read_csv(payload["scrna_matrix_path"], sep="\t", index_col=0)
        metadata = pd.read_csv(payload["scrna_metadata_path"], sep="\t")
        input_genes = Path(payload["input_gene_list_path"]).read_text().splitlines()
        output_genes = Path(payload["output_gene_list_path"]).read_text().splitlines()
        assert scrna.index.tolist() == ["Gene1", "Gene2"]
        assert metadata["cell_id"].tolist() == ["cell_1", "cell_2"]
        assert metadata["library_id"].tolist() == ["d1", "d2"]
        assert input_genes == ["Gene1", "Gene2"]
        assert output_genes == ["Gene1", "Gene2"]
        pd.DataFrame({"A": [0.2], "B": [0.8]}, index=["bulk_1"]).to_csv(payload["output_path"], sep="\t")
        Path(payload["ctgep_output_path"]).write_text("fake h5ad", encoding="utf-8")

    monkeypatch.setattr("foli_decon.tools.scrna_tools.run_command", fake_run_command)

    mixture = pd.DataFrame({"bulk_1": [10.0, 20.0]}, index=["Gene1", "Gene2"])
    scrna = pd.DataFrame(
        {
            "cell_1": [1.0, 2.0],
            "cell_2": [2.0, 1.0],
        },
        index=["Gene1", "Gene2"],
    )
    cell_types = pd.Series(["A", "B"], index=scrna.columns)
    batch_ids = pd.Series(["d1", "d2"], index=scrna.columns)

    result = run_deconvolution(
        tool="blue",
        mixture=mixture,
        scrna_counts=scrna,
        cell_types=cell_types,
        batch_ids=batch_ids,
        python_env="foli-decon-py311-cpu",
        blue_repo_path="/tmp/BLUE",
        n_training_samples=8,
        n_validation_samples=4,
        cells_per_sample=3,
        epochs=1,
        batch_size=2,
        extra_output_dir=str(tmp_path / "blue_extra"),
    )

    assert calls[0][:5] == ["mamba", "run", "-n", "foli-decon-py311-cpu", "python"]
    assert result.proportion.loc["bulk_1", "B"] == 0.8
    assert result.metadata["python_env"] == "foli-decon-py311-cpu"
    assert result.metadata["input_gene_count"] == 2
    assert result.output_paths["ctgep"] == str(tmp_path / "blue_extra" / "blue_predicted_ctGEP.h5ad")


@pytest.mark.parametrize(
    "tool_name,needle",
    [
        ("timer", "TIMER"),
        ("timer2", "TIMER2"),
        ("timer3", "TIMER3"),
        ("dsa", "DSA"),
    ],
)
def test_run_unsupported_requested_tool_raises_value_error_with_guidance(
    tool_name: str, needle: str
) -> None:
    """Requested but unavailable methods should fail with explicit environment guidance."""

    with pytest.raises(ValueError, match=needle):
        run_deconvolution(
            tool=tool_name,
            mixture=pd.DataFrame({"bulk_1": [1.0, 2.0]}, index=["Gene1", "Gene2"]),
            transform="raw",
        )
