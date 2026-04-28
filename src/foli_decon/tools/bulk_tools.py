"""Wrappers for bulk/signature deconvolution tools."""

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
from foli_decon.r_runner import run_r_tool


def _resolve_signature_matrix(
    signature: pd.DataFrame | None,
    scrna_counts: pd.DataFrame | None,
    cell_types: pd.Series | None,
    transform: str,
    gene_lengths: pd.Series | None,
) -> pd.DataFrame:
    """Resolve a signature matrix from either direct input or scRNA-seq aggregation."""

    if signature is not None:
        validate_expression_matrix(signature, "signature")
        return signature
    if scrna_counts is None or cell_types is None:
        raise ValueError("Provide either signature or both scrna_counts and cell_types")
    return build_signature_from_scrna(
        scrna_counts=scrna_counts,
        cell_types=cell_types,
        transform=transform,
        gene_lengths=gene_lengths,
    )


def run_xcell(
    mixture: pd.DataFrame,
    transform: str = "cpm",
    arrays: bool = False,
    expected_cell_types: list[str] | None = None,
) -> DeconvolutionResult:
    """Run xCell through the immunedeconv wrapper."""

    validate_expression_matrix(mixture, "mixture")
    mixture_norm = normalize_gene_names(mixture)
    mixture_transformed = apply_transform(mixture_norm, transform=transform)

    with tempfile.TemporaryDirectory(prefix="foli_xcell_") as workdir:
        workdir_path = Path(workdir)
        mixture_path = workdir_path / "mixture.tsv"
        write_expression_matrix(mixture_transformed, str(mixture_path))
        result, _ = run_r_tool(
            tool="xcell",
            payload={
                "mixture_path": str(mixture_path),
                "arrays": arrays,
                "expected_cell_types": expected_cell_types,
            },
        )

    return DeconvolutionResult(
        tool="xcell",
        proportions=result,
        metadata={
            "transform": transform,
            "arrays": arrays,
        },
    )


def run_xcell2(
    mixture: pd.DataFrame,
    xcell2_object_path: str,
    transform: str = "cpm",
    min_shared_genes: float = 0.2,
    raw_scores: bool = False,
    spillover: bool = True,
    spillover_alpha: float = 0.5,
) -> DeconvolutionResult:
    """Run xCell2 against a pre-trained xCell2 reference object."""

    validate_expression_matrix(mixture, "mixture")
    mixture_norm = normalize_gene_names(mixture)
    mixture_transformed = apply_transform(mixture_norm, transform=transform)
    xcell2_analysis_path = Path(__file__).resolve().parents[1] / "r_scripts" / "xCell2Analysis.R"

    with tempfile.TemporaryDirectory(prefix="foli_xcell2_") as workdir:
        workdir_path = Path(workdir)
        mixture_path = workdir_path / "mixture.tsv"
        write_expression_matrix(mixture_transformed, str(mixture_path))
        result, _ = run_r_tool(
            tool="xcell2",
            payload={
                "mixture_path": str(mixture_path),
                "xcell2_object_path": xcell2_object_path,
                "xcell2_analysis_path": str(xcell2_analysis_path),
                "min_shared_genes": min_shared_genes,
                "raw_scores": raw_scores,
                "spillover": spillover,
                "spillover_alpha": spillover_alpha,
            },
        )

    return DeconvolutionResult(
        tool="xcell2",
        proportions=result,
        metadata={
            "transform": transform,
            "min_shared_genes": min_shared_genes,
            "raw_scores": raw_scores,
            "spillover": spillover,
            "spillover_alpha": spillover_alpha,
        },
    )


def run_epic(
    mixture: pd.DataFrame,
    signature: pd.DataFrame | None = None,
    scrna_counts: pd.DataFrame | None = None,
    cell_types: pd.Series | None = None,
    mixture_transform: str = "cpm",
    signature_transform: str = "cpm",
    gene_lengths: pd.Series | None = None,
    tumor: bool = False,
    scale_mrna: bool = True,
) -> DeconvolutionResult:
    """Run EPIC with either built-in references or a custom signature."""

    validate_expression_matrix(mixture, "mixture")
    mixture_norm = normalize_gene_names(mixture)
    mixture_transformed = apply_transform(
        mixture_norm,
        transform=mixture_transform,
        gene_lengths=gene_lengths,
    )

    signature_matrix = None
    if signature is not None or scrna_counts is not None:
        signature_raw = _resolve_signature_matrix(
            signature=signature,
            scrna_counts=scrna_counts,
            cell_types=cell_types,
            transform=signature_transform,
            gene_lengths=gene_lengths,
        )
        signature_norm = normalize_gene_names(signature_raw)
        signature_transformed = apply_transform(
            signature_norm,
            transform=signature_transform,
            gene_lengths=gene_lengths,
        )
        mixture_transformed, signature_matrix = align_mixture_and_signature(
            mixture_transformed,
            signature_transformed,
        )

    with tempfile.TemporaryDirectory(prefix="foli_epic_") as workdir:
        workdir_path = Path(workdir)
        mixture_path = workdir_path / "mixture.tsv"
        write_expression_matrix(mixture_transformed, str(mixture_path))

        payload: dict[str, object] = {
            "mixture_path": str(mixture_path),
            "tumor": tumor,
            "scale_mrna": scale_mrna,
        }

        if signature_matrix is not None:
            signature_path = workdir_path / "signature.tsv"
            write_expression_matrix(signature_matrix, str(signature_path))
            payload["signature_path"] = str(signature_path)

        result, _ = run_r_tool(tool="epic", payload=payload)

    return DeconvolutionResult(
        tool="epic",
        proportions=result,
        metadata={
            "mixture_transform": mixture_transform,
            "signature_transform": signature_transform,
            "tumor": tumor,
            "scale_mrna": scale_mrna,
        },
    )


def run_dtangle(
    mixture: pd.DataFrame,
    signature: pd.DataFrame | None = None,
    scrna_counts: pd.DataFrame | None = None,
    cell_types: pd.Series | None = None,
    mixture_transform: str = "cpm_log1p",
    signature_transform: str = "cpm_log1p",
    gene_lengths: pd.Series | None = None,
    n_markers: int = 50,
    marker_method: str = "ratio",
) -> DeconvolutionResult:
    """Run dtangle with a custom signature matrix."""

    validate_expression_matrix(mixture, "mixture")
    signature_raw = _resolve_signature_matrix(
        signature=signature,
        scrna_counts=scrna_counts,
        cell_types=cell_types,
        transform=signature_transform,
        gene_lengths=gene_lengths,
    )

    mixture_norm = normalize_gene_names(mixture)
    signature_norm = normalize_gene_names(signature_raw)
    mixture_transformed = apply_transform(
        mixture_norm,
        transform=mixture_transform,
        gene_lengths=gene_lengths,
    )
    signature_transformed = apply_transform(
        signature_norm,
        transform=signature_transform,
        gene_lengths=gene_lengths,
    )
    mixture_aligned, signature_aligned = align_mixture_and_signature(
        mixture_transformed,
        signature_transformed,
    )

    with tempfile.TemporaryDirectory(prefix="foli_dtangle_") as workdir:
        workdir_path = Path(workdir)
        mixture_path = workdir_path / "mixture.tsv"
        signature_path = workdir_path / "signature.tsv"
        write_expression_matrix(mixture_aligned, str(mixture_path))
        write_expression_matrix(signature_aligned, str(signature_path))
        result, _ = run_r_tool(
            tool="dtangle",
            payload={
                "mixture_path": str(mixture_path),
                "signature_path": str(signature_path),
                "n_markers": n_markers,
                "marker_method": marker_method,
            },
        )

    return DeconvolutionResult(
        tool="dtangle",
        proportions=result,
        metadata={
            "mixture_transform": mixture_transform,
            "signature_transform": signature_transform,
            "n_markers": n_markers,
            "marker_method": marker_method,
        },
    )


def run_deconrnaseq(
    mixture: pd.DataFrame,
    signature: pd.DataFrame | None = None,
    scrna_counts: pd.DataFrame | None = None,
    cell_types: pd.Series | None = None,
    mixture_transform: str = "cpm",
    signature_transform: str = "cpm",
    gene_lengths: pd.Series | None = None,
    checksig: bool = False,
    use_scale: bool = True,
) -> DeconvolutionResult:
    """Run DeconRNASeq with a custom signature matrix."""

    validate_expression_matrix(mixture, "mixture")
    signature_raw = _resolve_signature_matrix(
        signature=signature,
        scrna_counts=scrna_counts,
        cell_types=cell_types,
        transform=signature_transform,
        gene_lengths=gene_lengths,
    )

    mixture_norm = normalize_gene_names(mixture)
    signature_norm = normalize_gene_names(signature_raw)
    mixture_transformed = apply_transform(
        mixture_norm,
        transform=mixture_transform,
        gene_lengths=gene_lengths,
    )
    signature_transformed = apply_transform(
        signature_norm,
        transform=signature_transform,
        gene_lengths=gene_lengths,
    )
    mixture_aligned, signature_aligned = align_mixture_and_signature(
        mixture_transformed,
        signature_transformed,
    )

    with tempfile.TemporaryDirectory(prefix="foli_deconrnaseq_") as workdir:
        workdir_path = Path(workdir)
        mixture_path = workdir_path / "mixture.tsv"
        signature_path = workdir_path / "signature.tsv"
        write_expression_matrix(mixture_aligned, str(mixture_path))
        write_expression_matrix(signature_aligned, str(signature_path))
        result, _ = run_r_tool(
            tool="deconrnaseq",
            payload={
                "mixture_path": str(mixture_path),
                "signature_path": str(signature_path),
                "checksig": checksig,
                "use_scale": use_scale,
            },
        )

    return DeconvolutionResult(
        tool="deconrnaseq",
        proportions=result,
        metadata={
            "mixture_transform": mixture_transform,
            "signature_transform": signature_transform,
            "checksig": checksig,
            "use_scale": use_scale,
        },
    )


def run_dwls(
    mixture: pd.DataFrame,
    signature: pd.DataFrame | None = None,
    scrna_counts: pd.DataFrame | None = None,
    cell_types: pd.Series | None = None,
    mixture_transform: str = "cpm",
    signature_transform: str = "cpm",
    gene_lengths: pd.Series | None = None,
    dwls_submethod: str = "DampenedWLS",
) -> DeconvolutionResult:
    """Run DWLS with either a direct signature or one built from scRNA-seq."""

    validate_expression_matrix(mixture, "mixture")
    signature_raw = _resolve_signature_matrix(
        signature=signature,
        scrna_counts=scrna_counts,
        cell_types=cell_types,
        transform=signature_transform,
        gene_lengths=gene_lengths,
    )

    mixture_norm = normalize_gene_names(mixture)
    signature_norm = normalize_gene_names(signature_raw)
    mixture_transformed = apply_transform(
        mixture_norm,
        transform=mixture_transform,
        gene_lengths=gene_lengths,
    )
    signature_transformed = apply_transform(
        signature_norm,
        transform=signature_transform,
        gene_lengths=gene_lengths,
    )
    mixture_aligned, signature_aligned = align_mixture_and_signature(
        mixture_transformed,
        signature_transformed,
    )

    with tempfile.TemporaryDirectory(prefix="foli_dwls_") as workdir:
        workdir_path = Path(workdir)
        mixture_path = workdir_path / "mixture.tsv"
        signature_path = workdir_path / "signature.tsv"
        write_expression_matrix(mixture_aligned, str(mixture_path))
        write_expression_matrix(signature_aligned, str(signature_path))
        result, _ = run_r_tool(
            tool="dwls",
            payload={
                "mixture_path": str(mixture_path),
                "signature_path": str(signature_path),
                "dwls_submethod": dwls_submethod,
            },
        )

    return DeconvolutionResult(
        tool="dwls",
        proportions=result,
        metadata={
            "mixture_transform": mixture_transform,
            "signature_transform": signature_transform,
            "dwls_submethod": dwls_submethod,
        },
    )
