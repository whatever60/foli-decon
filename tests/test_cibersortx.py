"""Unit tests for the CIBERSORTx Python wrapper."""

from pathlib import Path

import pandas as pd

from foli_decon.tools.cibersortx import run_cibersortx


def test_run_cibersortx_with_mocked_container(monkeypatch) -> None:
    """CIBERSORTx wrapper should return only cell-type columns from result output."""

    mixture = pd.DataFrame(
        {
            "sample_1": [10.0, 20.0, 30.0],
            "sample_2": [12.0, 22.0, 32.0],
        },
        index=["GeneA", "GeneB", "GeneC"],
    )
    signature = pd.DataFrame(
        {
            "CT_A": [50.0, 10.0, 5.0],
            "CT_B": [5.0, 45.0, 10.0],
        },
        index=["GeneA", "GeneB", "GeneC"],
    )

    captured: dict[str, list[str]] = {}

    def fake_run(command: list[str], check: bool, **kwargs) -> None:
        """Simulate a successful CIBERSORTx container run by writing expected output."""

        captured["command"] = command
        assert check is True
        out_mount = [value for value in command if value.endswith(":/src/outdir")][0]
        out_dir = Path(out_mount.replace(":/src/outdir", ""))

        simulated = pd.DataFrame(
            {
                "CT_A": [0.7, 0.4],
                "CT_B": [0.3, 0.6],
                "Correlation": [0.95, 0.90],
                "RMSE": [0.02, 0.03],
                "P-value": [0.001, 0.005],
            },
            index=["sample_1", "sample_2"],
        )
        simulated.to_csv(out_dir / "CIBERSORTx_Results.txt", sep="\t")

    monkeypatch.setattr("foli_decon.tools.cibersortx.subprocess.run", fake_run)

    result = run_cibersortx(
        mixture=mixture,
        signature=signature,
        username="user",
        token="token",
        mixture_transform="raw",
        signature_transform="raw",
        container_runtime="podman",
    )

    assert result.tool == "cibersortx"
    assert list(result.proportion.columns) == ["CT_A", "CT_B"]
    assert result.proportion.shape == (2, 2)
    assert captured["command"][0] == "podman"
