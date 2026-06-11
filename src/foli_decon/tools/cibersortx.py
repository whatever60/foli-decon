"""Wrapper for CIBERSORTx container execution."""

import subprocess
import tempfile
from pathlib import Path

import pandas as pd

from foli_decon.io import write_expression_matrix
from foli_decon.models import DeconvolutionResult
from foli_decon.preprocess import (
    align_mixture_and_signature,
    apply_transform,
    build_signature_from_scrna,
    normalize_gene_names,
    validate_expression_matrix,
)


def _redact_sensitive_text(text: str | None, username: str, token: str) -> str:
    """Redact CIBERSORTx credentials from process output."""

    if text is None:
        return ""
    redacted = text.replace(username, "<cibersortx_username>")
    redacted = redacted.replace(token, "<cibersortx_token>")
    return redacted


def _resolve_signature_matrix(
    signature: pd.DataFrame | None,
    scrna_counts: pd.DataFrame | None,
    cell_types: pd.Series | None,
    transform: str,
) -> pd.DataFrame:
    """Resolve a signature matrix from direct input or scRNA aggregation."""

    if signature is not None:
        validate_expression_matrix(signature, "signature")
        return signature
    if scrna_counts is None or cell_types is None:
        raise ValueError("Provide either signature or both scrna_counts and cell_types")
    return build_signature_from_scrna(
        scrna_counts=scrna_counts,
        cell_types=cell_types,
        transform=transform,
    )


def run_cibersortx(
    mixture: pd.DataFrame,
    username: str,
    token: str,
    signature: pd.DataFrame | None = None,
    scrna_counts: pd.DataFrame | None = None,
    cell_types: pd.Series | None = None,
    mixture_transform: str = "raw",
    signature_transform: str = "raw",
    container_runtime: str = "podman",
    image: str = "cibersortx/fractions",
    label: str | None = None,
    container_runtime_args: list[str] | None = None,
    timeout_seconds: int | None = None,
) -> DeconvolutionResult:
    """Run CIBERSORTx fractions mode inside a container runtime."""

    validate_expression_matrix(mixture, "mixture")
    signature_raw = _resolve_signature_matrix(
        signature=signature,
        scrna_counts=scrna_counts,
        cell_types=cell_types,
        transform=signature_transform,
    )

    mixture_norm = normalize_gene_names(mixture)
    signature_norm = normalize_gene_names(signature_raw)
    mixture_transformed = apply_transform(mixture_norm, transform=mixture_transform)
    signature_transformed = apply_transform(signature_norm, transform=signature_transform)
    mixture_aligned, signature_aligned = align_mixture_and_signature(
        mixture_transformed,
        signature_transformed,
    )

    with tempfile.TemporaryDirectory(prefix="foli_cibersortx_") as workdir:
        workdir_path = Path(workdir)
        output_path = workdir_path / "output"
        output_path.mkdir(parents=True, exist_ok=True)

        mixture_path = workdir_path / "mixture.tsv"
        signature_path = workdir_path / "sigmatrix.tsv"

        write_expression_matrix(mixture_aligned, str(mixture_path))
        write_expression_matrix(signature_aligned, str(signature_path))

        command: list[str] = [
            container_runtime,
            "run",
            "--rm",
        ]
        runtime_args = container_runtime_args
        if runtime_args is None and container_runtime == "podman":
            runtime_args = ["--runtime=runc", "--network=slirp4netns"]
        if runtime_args is not None:
            command.extend(runtime_args)
        command.extend(
            [
                "-v",
                f"{workdir}:/src/data",
                "-v",
                f"{output_path}:/src/outdir",
                image,
                "--username",
                username,
                "--token",
                token,
                "--mixture",
                "mixture.tsv",
                "--sigmatrix",
                "sigmatrix.tsv",
            ]
        )
        if label is not None:
            command.extend(["--label", label])

        try:
            subprocess.run(
                command,
                check=True,
                timeout=timeout_seconds,
                capture_output=True,
                text=True,
            )
        except subprocess.TimeoutExpired as error:
            raise RuntimeError(f"CIBERSORTx fractions timed out after {timeout_seconds} seconds") from error
        except subprocess.CalledProcessError as error:
            stdout = _redact_sensitive_text(error.stdout, username=username, token=token)
            stderr = _redact_sensitive_text(error.stderr, username=username, token=token)
            raise RuntimeError(
                f"CIBERSORTx fractions failed with exit code {error.returncode}.\n"
                f"stdout:\n{stdout[-1200:]}\n"
                f"stderr:\n{stderr[-1200:]}"
            ) from error

        if label is None:
            result_file = output_path / "CIBERSORTx_Results.txt"
        else:
            result_file = output_path / f"CIBERSORTx_{label}_Results.txt"

        result = pd.read_csv(result_file, sep="\t", index_col=0)
        extra_columns = ["Correlation", "RMSE", "P-value"]
        kept_columns = [column for column in result.columns if column not in extra_columns]
        result = result[kept_columns]

    return DeconvolutionResult(
        tool="cibersortx",
        proportion=result,
        metadata={
            "mixture_transform": mixture_transform,
            "signature_transform": signature_transform,
            "container_runtime": container_runtime,
            "image": image,
            "label_used": label is not None,
            "container_runtime_args": runtime_args,
            "timeout_seconds": 0 if timeout_seconds is None else timeout_seconds,
        },
    )
