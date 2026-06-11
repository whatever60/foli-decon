"""Wrappers for bulk/signature deconvolution tools."""

import json
import tempfile
from pathlib import Path

import numpy as np
import pandas as pd

from foli_decon.io import read_sample_by_celltype, write_expression_matrix, write_vector
from foli_decon.models import DeconvolutionResult, FeatureSelectionResult, ReferenceTrainingResult
from foli_decon.preprocess import (
    align_mixture_and_signature,
    apply_transform,
    build_signature_from_scrna,
    intersect_gene_sets,
    normalize_gene_names,
    validate_expression_matrix,
)
from foli_decon.python_runner import run_command
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


def _prepare_scrna_feature_selection(
    scrna_counts: pd.DataFrame,
    cell_types: pd.Series,
    scrna_transform: str,
    gene_lengths: pd.Series | None,
    nfeatures: int,
) -> tuple[pd.DataFrame, pd.Series]:
    """Validate and transform scRNA data for label-aware feature selectors."""

    validate_expression_matrix(scrna_counts, "scrna_counts")
    if len(cell_types) != scrna_counts.shape[1]:
        raise ValueError("cell_types length must match number of single-cell columns")
    scrna_norm = normalize_gene_names(scrna_counts)
    scrna_transformed = apply_transform(
        scrna_norm,
        transform=scrna_transform,
        gene_lengths=gene_lengths,
    )
    if nfeatures >= scrna_transformed.shape[0]:
        raise ValueError("nfeatures must be smaller than the number of genes for this native feature selector")
    labels = pd.Series(cell_types.to_numpy(), index=scrna_transformed.columns).astype(str)
    return scrna_transformed, labels


def run_xcell(
    mixture: pd.DataFrame,
    transform: str = "cpm",
    arrays: bool = False,
    expected_cell_types: list[str] | None = None,
    min_genes: int = 5000,
) -> DeconvolutionResult:
    """Run xCell through its native R API."""

    validate_expression_matrix(mixture, "mixture")
    mixture_norm = normalize_gene_names(mixture)
    mixture_transformed = apply_transform(mixture_norm, transform=transform)
    if mixture_transformed.shape[0] < min_genes:
        raise ValueError(
            "xCell requires at least "
            f"{min_genes} input genes before its internal gene filtering; got {mixture_transformed.shape[0]}. "
            "Built-in xCell is not suitable for targeted-panel inputs this small."
        )

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
        score=result,
        metadata={
            "transform": transform,
            "arrays": arrays,
            "min_genes": min_genes,
        },
    )


def run_quantiseq(
    mixture: pd.DataFrame,
    transform: str = "cpm",
    tumor: bool = False,
    arrays: bool = False,
    scale_mrna: bool = True,
) -> DeconvolutionResult:
    """Run quanTIseq through its native R API."""

    validate_expression_matrix(mixture, "mixture")
    mixture_norm = normalize_gene_names(mixture)
    mixture_transformed = apply_transform(mixture_norm, transform=transform)

    with tempfile.TemporaryDirectory(prefix="foli_quantiseq_") as workdir:
        workdir_path = Path(workdir)
        mixture_path = workdir_path / "mixture.tsv"
        write_expression_matrix(mixture_transformed, str(mixture_path))
        result, _ = run_r_tool(
            tool="quantiseq",
            payload={
                "mixture_path": str(mixture_path),
                "tumor": tumor,
                "arrays": arrays,
                "scale_mrna": scale_mrna,
            },
        )

    return DeconvolutionResult(
        tool="quantiseq",
        proportion=result,
        metadata={
            "transform": transform,
            "tumor": tumor,
            "arrays": arrays,
            "scale_mrna": scale_mrna,
        },
    )


def run_mcp_counter(
    mixture: pd.DataFrame,
    transform: str = "cpm",
    feature_types: str = "HUGO_symbols",
    probesets_path: str | None = None,
    genes_path: str | None = None,
) -> DeconvolutionResult:
    """Run MCP-counter through its native R API."""

    validate_expression_matrix(mixture, "mixture")
    mixture_norm = normalize_gene_names(mixture)
    mixture_transformed = apply_transform(mixture_norm, transform=transform)

    with tempfile.TemporaryDirectory(prefix="foli_mcp_counter_") as workdir:
        workdir_path = Path(workdir)
        mixture_path = workdir_path / "mixture.tsv"
        write_expression_matrix(mixture_transformed, str(mixture_path))
        result, _ = run_r_tool(
            tool="mcp_counter",
            payload={
                "mixture_path": str(mixture_path),
                "feature_types": feature_types,
                "probesets_path": "" if probesets_path is None else probesets_path,
                "genes_path": "" if genes_path is None else genes_path,
            },
        )

    return DeconvolutionResult(
        tool="mcp_counter",
        score=result,
        metadata={
            "transform": transform,
            "feature_types": feature_types,
            "probesets_path": "" if probesets_path is None else probesets_path,
            "genes_path": "" if genes_path is None else genes_path,
        },
    )


def run_psea(
    mixture: pd.DataFrame,
    marker_sets: dict[str, list[str] | list[list[str]]],
    transform: str = "cpm",
    sample_subset: list[str] | None = None,
    target_mean: float = 1.0,
) -> DeconvolutionResult:
    """Run PSEA marker reference-signal scoring."""

    validate_expression_matrix(mixture, "mixture")
    mixture_norm = normalize_gene_names(mixture)
    mixture_transformed = apply_transform(mixture_norm, transform=transform)

    marker_rows = []
    for cell_type, marker_groups in marker_sets.items():
        for group_index, marker_group in enumerate(marker_groups):
            if isinstance(marker_group, str):
                genes = [marker_group]
            else:
                genes = list(marker_group)
            for gene in genes:
                marker_rows.append(
                    {
                        "cell_type": cell_type,
                        "marker_group": f"{cell_type}_{group_index}",
                        "gene": str(gene).replace(".", "-"),
                    }
                )
    marker_table = pd.DataFrame(marker_rows)

    with tempfile.TemporaryDirectory(prefix="foli_psea_") as workdir:
        workdir_path = Path(workdir)
        mixture_path = workdir_path / "mixture.tsv"
        marker_sets_path = workdir_path / "marker_sets.tsv"
        write_expression_matrix(mixture_transformed, str(mixture_path))
        marker_table.to_csv(marker_sets_path, sep="\t", index=False)
        result, _ = run_r_tool(
            tool="psea",
            payload={
                "mixture_path": str(mixture_path),
                "marker_sets_path": str(marker_sets_path),
                "sample_subset": sample_subset,
                "target_mean": target_mean,
            },
        )

    return DeconvolutionResult(
        tool="psea",
        score=result,
        metadata={
            "transform": transform,
            "n_marker_rows": int(marker_table.shape[0]),
            "n_cell_types": int(marker_table["cell_type"].nunique()),
            "sample_subset_count": 0 if sample_subset is None else len(sample_subset),
            "target_mean": target_mean,
            "output_semantics": "PSEA marker reference signal; not a compositional fraction",
        },
    )


def run_cellcode(
    mixture: pd.DataFrame,
    signature: pd.DataFrame | None = None,
    scrna_counts: pd.DataFrame | None = None,
    cell_types: pd.Series | None = None,
    groups: pd.Series | None = None,
    mixture_transform: str = "cpm_log1p",
    signature_transform: str = "cpm_log1p",
    gene_lengths: pd.Series | None = None,
    tag_cutoff: float = 2.0,
    max_markers: int | None = None,
    ref_mean: bool = False,
    method: str = "mixed",
    mix_par: float = 0.3,
    r_env: str | None = None,
    r_output_path: str | None = None,
    extra_output_dir: str | None = None,
) -> DeconvolutionResult:
    """Run CellCODE surrogate proportion variable scoring."""

    validate_expression_matrix(mixture, "mixture")
    if groups is None:
        raise ValueError("groups is required for CellCODE because the native API models sample groups")
    if len(groups) != mixture.shape[1]:
        raise ValueError("groups length must match number of mixture samples")
    if pd.Series(groups.to_numpy()).nunique() < 2:
        raise ValueError("CellCODE groups must contain at least two levels")

    if signature is None:
        signature_raw = _resolve_signature_matrix(
            signature=signature,
            scrna_counts=scrna_counts,
            cell_types=cell_types,
            transform=signature_transform,
            gene_lengths=gene_lengths,
        )
        signature_needs_transform = False
    else:
        signature_raw = _resolve_signature_matrix(
            signature=signature,
            scrna_counts=scrna_counts,
            cell_types=cell_types,
            transform=signature_transform,
            gene_lengths=gene_lengths,
        )
        signature_needs_transform = True
    mixture_norm = normalize_gene_names(mixture)
    signature_norm = normalize_gene_names(signature_raw)
    mixture_transformed = apply_transform(mixture_norm, transform=mixture_transform, gene_lengths=gene_lengths)
    signature_transformed = (
        apply_transform(signature_norm, transform=signature_transform, gene_lengths=gene_lengths)
        if signature_needs_transform
        else signature_norm
    )
    mixture_aligned, signature_aligned = intersect_gene_sets([mixture_transformed, signature_transformed])

    rscript_command = None
    if r_env is not None:
        rscript_command = ["mamba", "run", "-n", r_env, "Rscript"]

    with tempfile.TemporaryDirectory(prefix="foli_cellcode_") as workdir:
        workdir_path = Path(workdir)
        mixture_path = workdir_path / "mixture.tsv"
        signature_path = workdir_path / "signature.tsv"
        groups_path = workdir_path / "groups.tsv"
        write_expression_matrix(mixture_aligned, str(mixture_path))
        write_expression_matrix(signature_aligned, str(signature_path))
        write_vector(pd.Series(groups.to_numpy()), str(groups_path), "group")

        result, _ = run_r_tool(
            tool="cellcode",
            payload={
                "mixture_path": str(mixture_path),
                "signature_path": str(signature_path),
                "groups_path": str(groups_path),
                "tag_cutoff": tag_cutoff,
                "max_markers": max_markers,
                "ref_mean": ref_mean,
                "method": method,
                "mix_par": mix_par,
                "extra_output_dir": extra_output_dir,
            },
            rscript_command=rscript_command,
            output_path=r_output_path,
        )

    output_paths = {}
    if r_output_path is not None:
        output_paths["score"] = r_output_path
    if extra_output_dir is not None:
        output_paths["extra_output_dir"] = extra_output_dir

    return DeconvolutionResult(
        tool="cellcode",
        score=result,
        metadata={
            "mixture_transform": mixture_transform,
            "signature_transform": signature_transform,
            "tag_cutoff": tag_cutoff,
            "max_markers": 0 if max_markers is None else max_markers,
            "ref_mean": ref_mean,
            "method": method,
            "mix_par": mix_par,
            "r_env": "" if r_env is None else r_env,
            "r_output_path": "" if r_output_path is None else r_output_path,
            "extra_output_dir": "" if extra_output_dir is None else extra_output_dir,
            "output_semantics": "CellCODE surrogate proportion variables; scores, not compositional fractions",
        },
        output_paths=output_paths,
    )


def select_features_autogenes(
    signature: pd.DataFrame | None = None,
    scrna_counts: pd.DataFrame | None = None,
    cell_types: pd.Series | None = None,
    signature_transform: str = "cpm",
    gene_lengths: pd.Series | None = None,
    python_env: str | None = None,
    ngen: int = 50,
    nfeatures: int = 500,
    mode: str = "fixed",
    seed: int = 0,
    selection_index: int = 0,
    verbose: bool = False,
    selected_features_path: str | None = None,
    log_path: str | None = None,
) -> FeatureSelectionResult:
    """Run AutoGeneS feature selection through a Python sidecar."""

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

    command_base = ["python"] if python_env is None else ["mamba", "run", "-n", python_env, "python"]
    runner_path = Path(__file__).resolve().parents[1] / "deep_scripts" / "select_autogenes.py"

    with tempfile.TemporaryDirectory(prefix="foli_autogenes_select_") as workdir:
        workdir_path = Path(workdir)
        signature_path = workdir_path / "signature.tsv"
        payload_path = workdir_path / "payload.json"
        resolved_selected_features_path = (
            workdir_path / "autogenes_selected_features.tsv"
            if selected_features_path is None
            else Path(selected_features_path)
        )

        write_expression_matrix(signature_transformed, str(signature_path))
        payload = {
            "signature_path": str(signature_path),
            "selected_features_path": str(resolved_selected_features_path),
            "ngen": ngen,
            "nfeatures": nfeatures,
            "mode": mode,
            "seed": seed,
            "selection_index": selection_index,
            "verbose": verbose,
        }
        payload_path.write_text(json.dumps(payload), encoding="utf-8")
        run_command(command_base + [str(runner_path), str(payload_path)], log_path=log_path)
        selected_features = pd.read_csv(resolved_selected_features_path, sep="\t")

    output_paths = {}
    if selected_features_path is not None:
        output_paths["selected_features"] = selected_features_path

    return FeatureSelectionResult(
        tool="autogenes",
        selected_genes=pd.Index(selected_features["gene"].astype(str)),
        ranking=None,
        metadata={
            "signature_transform": signature_transform,
            "python_env": "" if python_env is None else python_env,
            "ngen": ngen,
            "nfeatures": nfeatures,
            "mode": mode,
            "seed": seed,
            "selection_index": selection_index,
            "log_path": "" if log_path is None else log_path,
        },
        output_paths=output_paths,
    )


def select_features_scgenefit(
    scrna_counts: pd.DataFrame,
    cell_types: pd.Series,
    scrna_transform: str = "cpm_log1p",
    gene_lengths: pd.Series | None = None,
    python_env: str | None = None,
    nfeatures: int = 50,
    method: str = "centers",
    epsilon: float = 1.0,
    sampling_rate: float = 1.0,
    n_neighbors: int = 3,
    max_constraints: int = 1000,
    redundancy: float = 0.01,
    seed: int = 0,
    verbose: bool = False,
    selected_features_path: str | None = None,
    ranking_path: str | None = None,
    log_path: str | None = None,
) -> FeatureSelectionResult:
    """Run scGeneFit marker selection through a Python sidecar."""

    scrna_transformed, labels = _prepare_scrna_feature_selection(
        scrna_counts=scrna_counts,
        cell_types=cell_types,
        scrna_transform=scrna_transform,
        gene_lengths=gene_lengths,
        nfeatures=nfeatures,
    )

    command_base = ["python"] if python_env is None else ["mamba", "run", "-n", python_env, "python"]
    runner_path = Path(__file__).resolve().parents[1] / "deep_scripts" / "select_scgenefit.py"

    with tempfile.TemporaryDirectory(prefix="foli_scgenefit_select_") as workdir:
        workdir_path = Path(workdir)
        scrna_path = workdir_path / "scrna.tsv"
        labels_path = workdir_path / "cell_types.tsv"
        payload_path = workdir_path / "payload.json"
        resolved_selected_features_path = (
            workdir_path / "scgenefit_selected_features.tsv"
            if selected_features_path is None
            else Path(selected_features_path)
        )
        resolved_ranking_path = workdir_path / "scgenefit_ranking.tsv" if ranking_path is None else Path(ranking_path)

        write_expression_matrix(scrna_transformed, str(scrna_path))
        write_vector(labels, str(labels_path), "cell_type")
        payload = {
            "scrna_path": str(scrna_path),
            "labels_path": str(labels_path),
            "selected_features_path": str(resolved_selected_features_path),
            "ranking_path": str(resolved_ranking_path),
            "nfeatures": nfeatures,
            "method": method,
            "epsilon": epsilon,
            "sampling_rate": sampling_rate,
            "n_neighbors": n_neighbors,
            "max_constraints": max_constraints,
            "redundancy": redundancy,
            "seed": seed,
            "verbose": verbose,
        }
        payload_path.write_text(json.dumps(payload), encoding="utf-8")
        run_command(command_base + [str(runner_path), str(payload_path)], log_path=log_path)
        selected_features = pd.read_csv(resolved_selected_features_path, sep="\t")
        ranking = pd.read_csv(resolved_ranking_path, sep="\t")

    output_paths = {}
    if selected_features_path is not None:
        output_paths["selected_features"] = selected_features_path
    if ranking_path is not None:
        output_paths["ranking"] = ranking_path

    return FeatureSelectionResult(
        tool="scgenefit",
        selected_genes=pd.Index(selected_features["gene"].astype(str)),
        ranking=ranking,
        metadata={
            "scrna_transform": scrna_transform,
            "python_env": "" if python_env is None else python_env,
            "nfeatures": nfeatures,
            "method": method,
            "epsilon": epsilon,
            "sampling_rate": sampling_rate,
            "n_neighbors": n_neighbors,
            "max_constraints": max_constraints,
            "redundancy": redundancy,
            "seed": seed,
            "log_path": "" if log_path is None else log_path,
            "native_input": "dense cells x genes matrix plus labels",
        },
        output_paths=output_paths,
    )


def select_features_markermap(
    scrna_counts: pd.DataFrame,
    cell_types: pd.Series,
    scrna_transform: str = "cpm_log1p",
    gene_lengths: pd.Series | None = None,
    python_env: str | None = None,
    nfeatures: int = 50,
    hidden_layer_size: int | None = None,
    z_size: int = 16,
    batch_size: int = 64,
    loss_tradeoff: float = 0.0,
    train_fraction: float = 0.8,
    min_epochs: int = 2,
    max_epochs: int = 10,
    auto_lr: bool = False,
    early_stopping_patience: int = 3,
    seed: int = 0,
    verbose: bool = False,
    selected_features_path: str | None = None,
    ranking_path: str | None = None,
    log_path: str | None = None,
) -> FeatureSelectionResult:
    """Run MarkerMap supervised marker selection through a Python sidecar."""

    scrna_transformed, labels = _prepare_scrna_feature_selection(
        scrna_counts=scrna_counts,
        cell_types=cell_types,
        scrna_transform=scrna_transform,
        gene_lengths=gene_lengths,
        nfeatures=nfeatures,
    )
    resolved_hidden_layer_size = (
        max(4, min(64, int(np.ceil(scrna_transformed.shape[0] * 0.1))))
        if hidden_layer_size is None
        else hidden_layer_size
    )

    command_base = ["python"] if python_env is None else ["mamba", "run", "-n", python_env, "python"]
    runner_path = Path(__file__).resolve().parents[1] / "deep_scripts" / "select_markermap.py"

    with tempfile.TemporaryDirectory(prefix="foli_markermap_select_") as workdir:
        workdir_path = Path(workdir)
        scrna_path = workdir_path / "scrna.tsv"
        labels_path = workdir_path / "cell_types.tsv"
        payload_path = workdir_path / "payload.json"
        resolved_selected_features_path = (
            workdir_path / "markermap_selected_features.tsv"
            if selected_features_path is None
            else Path(selected_features_path)
        )
        resolved_ranking_path = workdir_path / "markermap_ranking.tsv" if ranking_path is None else Path(ranking_path)

        write_expression_matrix(scrna_transformed, str(scrna_path))
        write_vector(labels, str(labels_path), "cell_type")
        payload = {
            "scrna_path": str(scrna_path),
            "labels_path": str(labels_path),
            "selected_features_path": str(resolved_selected_features_path),
            "ranking_path": str(resolved_ranking_path),
            "nfeatures": nfeatures,
            "hidden_layer_size": resolved_hidden_layer_size,
            "z_size": z_size,
            "batch_size": batch_size,
            "loss_tradeoff": loss_tradeoff,
            "train_fraction": train_fraction,
            "min_epochs": min_epochs,
            "max_epochs": max_epochs,
            "auto_lr": auto_lr,
            "early_stopping_patience": early_stopping_patience,
            "seed": seed,
            "verbose": verbose,
        }
        payload_path.write_text(json.dumps(payload), encoding="utf-8")
        run_command(command_base + [str(runner_path), str(payload_path)], log_path=log_path)
        selected_features = pd.read_csv(resolved_selected_features_path, sep="\t")
        ranking = pd.read_csv(resolved_ranking_path, sep="\t")

    output_paths = {}
    if selected_features_path is not None:
        output_paths["selected_features"] = selected_features_path
    if ranking_path is not None:
        output_paths["ranking"] = ranking_path

    return FeatureSelectionResult(
        tool="markermap",
        selected_genes=pd.Index(selected_features["gene"].astype(str)),
        ranking=ranking,
        metadata={
            "scrna_transform": scrna_transform,
            "python_env": "" if python_env is None else python_env,
            "nfeatures": nfeatures,
            "hidden_layer_size": resolved_hidden_layer_size,
            "z_size": z_size,
            "batch_size": batch_size,
            "loss_tradeoff": loss_tradeoff,
            "train_fraction": train_fraction,
            "min_epochs": min_epochs,
            "max_epochs": max_epochs,
            "auto_lr": auto_lr,
            "early_stopping_patience": early_stopping_patience,
            "seed": seed,
            "log_path": "" if log_path is None else log_path,
            "native_input": "AnnData cells x genes matrix with obs['cell_type'] labels",
        },
        output_paths=output_paths,
    )


def run_autogenes(
    mixture: pd.DataFrame,
    signature: pd.DataFrame | None = None,
    scrna_counts: pd.DataFrame | None = None,
    cell_types: pd.Series | None = None,
    mixture_transform: str = "cpm",
    signature_transform: str = "cpm",
    gene_lengths: pd.Series | None = None,
    python_env: str | None = None,
    ngen: int = 50,
    nfeatures: int = 500,
    mode: str = "fixed",
    seed: int = 0,
    selection_index: int = 0,
    model: str = "nusvr",
    clip_negative: bool = True,
    normalize_coefficients: bool = True,
    verbose: bool = False,
    output_path: str | None = None,
    selected_features_path: str | None = None,
    log_path: str | None = None,
) -> DeconvolutionResult:
    """Run AutoGeneS feature selection and deconvolution through a Python sidecar."""

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

    command_base = ["python"] if python_env is None else ["mamba", "run", "-n", python_env, "python"]
    runner_path = Path(__file__).resolve().parents[1] / "deep_scripts" / "run_autogenes.py"

    with tempfile.TemporaryDirectory(prefix="foli_autogenes_") as workdir:
        workdir_path = Path(workdir)
        mixture_path = workdir_path / "mixture.tsv"
        signature_path = workdir_path / "signature.tsv"
        payload_path = workdir_path / "payload.json"
        resolved_output_path = workdir_path / "autogenes_proportions.tsv" if output_path is None else Path(output_path)
        resolved_selected_features_path = (
            workdir_path / "autogenes_selected_features.tsv"
            if selected_features_path is None
            else Path(selected_features_path)
        )

        write_expression_matrix(mixture_aligned, str(mixture_path))
        write_expression_matrix(signature_aligned, str(signature_path))
        payload = {
            "mixture_path": str(mixture_path),
            "signature_path": str(signature_path),
            "output_path": str(resolved_output_path),
            "selected_features_path": str(resolved_selected_features_path),
            "ngen": ngen,
            "nfeatures": nfeatures,
            "mode": mode,
            "seed": seed,
            "selection_index": selection_index,
            "model": model,
            "clip_negative": clip_negative,
            "normalize_coefficients": normalize_coefficients,
            "verbose": verbose,
        }
        payload_path.write_text(json.dumps(payload), encoding="utf-8")
        run_command(command_base + [str(runner_path), str(payload_path)], log_path=log_path)
        result = read_sample_by_celltype(str(resolved_output_path))
        selected_features = pd.read_csv(resolved_selected_features_path, sep="\t")

    output_paths = {}
    if output_path is not None:
        output_paths["proportion"] = output_path
    if selected_features_path is not None:
        output_paths["selected_features"] = selected_features_path

    return DeconvolutionResult(
        tool="autogenes",
        proportion=result,
        selected_features=selected_features,
        metadata={
            "mixture_transform": mixture_transform,
            "signature_transform": signature_transform,
            "python_env": "" if python_env is None else python_env,
            "ngen": ngen,
            "nfeatures": nfeatures,
            "mode": mode,
            "seed": seed,
            "selection_index": selection_index,
            "model": model,
            "clip_negative": clip_negative,
            "normalize_coefficients": normalize_coefficients,
            "log_path": "" if log_path is None else log_path,
        },
        output_paths=output_paths,
    )


def run_blade(
    mixture: pd.DataFrame,
    signature: pd.DataFrame | None = None,
    signature_sd: pd.DataFrame | None = None,
    scrna_counts: pd.DataFrame | None = None,
    cell_types: pd.Series | None = None,
    mixture_transform: str = "raw",
    signature_transform: str = "cpm",
    gene_lengths: pd.Series | None = None,
    signature_is_log: bool = False,
    signature_sd_floor: float = 0.01,
    python_env: str | None = None,
    alphas: list[float] | None = None,
    alpha0s: list[float] | None = None,
    kappa0s: list[float] | None = None,
    sigma_ys: list[float] | None = None,
    n_rep: int = 3,
    n_jobs: int = 1,
    n_rep_final: int = 3,
    feature_selection_fraction: float = 0.0,
    seed: int = 0,
    output_path: str | None = None,
    extra_output_dir: str | None = None,
    log_path: str | None = None,
) -> DeconvolutionResult:
    """Run BLADE through the Python sidecar package."""

    validate_expression_matrix(mixture, "mixture")
    mixture_norm = normalize_gene_names(mixture)
    mixture_transformed = apply_transform(
        mixture_norm,
        transform=mixture_transform,
        gene_lengths=gene_lengths,
    )

    reference_mode = "signature"
    if scrna_counts is not None:
        if cell_types is None:
            raise ValueError("cell_types are required when scrna_counts is provided")
        validate_expression_matrix(scrna_counts, "scrna_counts")
        if len(cell_types) != scrna_counts.shape[1]:
            raise ValueError("cell_types length must match number of scRNA columns")
        scrna_norm = normalize_gene_names(scrna_counts)
        scrna_transformed = apply_transform(
            scrna_norm,
            transform=signature_transform,
            gene_lengths=gene_lengths,
        )
        mixture_aligned, scrna_aligned = intersect_gene_sets([mixture_transformed, scrna_transformed])
        aligned_cell_types = pd.Series(cell_types.to_numpy(), index=scrna_norm.columns).loc[scrna_aligned.columns]
        scrna_log = np.log1p(scrna_aligned)
        signature_mean = scrna_log.T.groupby(aligned_cell_types).mean().T
        signature_sd_matrix = scrna_log.T.groupby(aligned_cell_types).std().T
        signature_sd_matrix = signature_sd_matrix.fillna(signature_sd_floor).clip(lower=signature_sd_floor)
        reference_mode = "scrna"
    else:
        if signature is None:
            raise ValueError("Provide either signature or both scrna_counts and cell_types")
        validate_expression_matrix(signature, "signature")
        signature_norm = normalize_gene_names(signature)
        if signature_is_log:
            signature_mean_raw = signature_norm
        else:
            signature_transformed = apply_transform(
                signature_norm,
                transform=signature_transform,
                gene_lengths=gene_lengths,
            )
            signature_mean_raw = np.log1p(signature_transformed)
        if signature_sd is None:
            signature_sd_norm = pd.DataFrame(
                signature_sd_floor,
                index=signature_mean_raw.index,
                columns=signature_mean_raw.columns,
            )
        else:
            validate_expression_matrix(signature_sd, "signature_sd")
            signature_sd_norm = normalize_gene_names(signature_sd)
        mixture_aligned, signature_mean, signature_sd_matrix = intersect_gene_sets(
            [mixture_transformed, signature_mean_raw, signature_sd_norm]
        )
        signature_sd_matrix = signature_sd_matrix.fillna(signature_sd_floor).clip(lower=signature_sd_floor)

    resolved_alphas = [1.0, 10.0] if alphas is None else alphas
    resolved_alpha0s = [0.1, 1.0, 5.0] if alpha0s is None else alpha0s
    resolved_kappa0s = [1.0, 0.5, 0.1] if kappa0s is None else kappa0s
    resolved_sigma_ys = [1.0, 0.3, 0.5] if sigma_ys is None else sigma_ys
    command_base = ["python"] if python_env is None else ["mamba", "run", "-n", python_env, "python"]
    runner_path = Path(__file__).resolve().parents[1] / "deep_scripts" / "run_blade.py"

    with tempfile.TemporaryDirectory(prefix="foli_blade_") as workdir:
        workdir_path = Path(workdir)
        mixture_path = workdir_path / "mixture.tsv"
        signature_mean_path = workdir_path / "signature_mean.tsv"
        signature_sd_path = workdir_path / "signature_sd.tsv"
        payload_path = workdir_path / "payload.json"
        resolved_output_path = workdir_path / "blade_proportions.tsv" if output_path is None else Path(output_path)

        write_expression_matrix(mixture_aligned, str(mixture_path))
        write_expression_matrix(signature_mean, str(signature_mean_path))
        write_expression_matrix(signature_sd_matrix, str(signature_sd_path))
        payload = {
            "mixture_path": str(mixture_path),
            "signature_mean_path": str(signature_mean_path),
            "signature_sd_path": str(signature_sd_path),
            "output_path": str(resolved_output_path),
            "extra_output_dir": "" if extra_output_dir is None else extra_output_dir,
            "alphas": resolved_alphas,
            "alpha0s": resolved_alpha0s,
            "kappa0s": resolved_kappa0s,
            "sigma_ys": resolved_sigma_ys,
            "n_rep": n_rep,
            "n_jobs": n_jobs,
            "n_rep_final": n_rep_final,
            "feature_selection_fraction": feature_selection_fraction,
            "seed": seed,
        }
        payload_path.write_text(json.dumps(payload), encoding="utf-8")
        run_command(command_base + [str(runner_path), str(payload_path)], log_path=log_path)
        result = read_sample_by_celltype(str(resolved_output_path))

    output_paths = {}
    if output_path is not None:
        output_paths["proportion"] = output_path
    if extra_output_dir is not None:
        output_paths["extra_output_dir"] = extra_output_dir

    return DeconvolutionResult(
        tool="blade",
        proportion=result,
        metadata={
            "mixture_transform": mixture_transform,
            "signature_transform": signature_transform,
            "reference_mode": reference_mode,
            "signature_is_log": signature_is_log,
            "signature_sd_floor": signature_sd_floor,
            "python_env": "" if python_env is None else python_env,
            "n_rep": n_rep,
            "n_jobs": n_jobs,
            "n_rep_final": n_rep_final,
            "feature_selection_fraction": feature_selection_fraction,
            "seed": seed,
            "log_path": "" if log_path is None else log_path,
        },
        output_paths=output_paths,
    )


def _default_abis_resource_path(technology: str, resource_name: str) -> Path:
    """Return the default local ABIS resource path for the source checkout."""

    resource_dir = Path(__file__).resolve().parents[3] / "resources" / "abis"
    technology_dirs = {
        "rnaseq": resource_dir / "rnaseq",
        "microarray": resource_dir / "microarray",
    }
    return technology_dirs[technology] / resource_name


def run_abis(
    mixture: pd.DataFrame,
    transform: str = "cpm",
    technology: str = "rnaseq",
    signature: pd.DataFrame | None = None,
    signature_path: str | None = None,
    target_path: str | None = None,
    gene_lengths: pd.Series | None = None,
    rescale_percent: bool = True,
    r_env: str | None = None,
    output_path: str | None = None,
    log_path: str | None = None,
) -> DeconvolutionResult:
    """Run ABIS using the local Shiny-app source logic and signature matrices."""

    validate_expression_matrix(mixture, "mixture")
    technology_normalized = technology.lower()
    if technology_normalized not in ("rnaseq", "microarray"):
        raise ValueError("technology must be 'rnaseq' or 'microarray'")
    mixture_norm = normalize_gene_names(mixture)
    mixture_transformed = apply_transform(mixture_norm, transform=transform, gene_lengths=gene_lengths)

    rscript_command = None
    if r_env is not None:
        rscript_command = ["mamba", "run", "-n", r_env, "Rscript"]

    with tempfile.TemporaryDirectory(prefix="foli_abis_") as workdir:
        workdir_path = Path(workdir)
        mixture_path = workdir_path / "mixture.tsv"
        payload_signature_path = signature_path
        payload_target_path = target_path
        write_expression_matrix(mixture_transformed, str(mixture_path))

        if signature is not None:
            validate_expression_matrix(signature, "signature")
            signature_norm = normalize_gene_names(signature)
            payload_signature_path = str(workdir_path / "abis_signature.tsv")
            write_expression_matrix(signature_norm, payload_signature_path)
        if payload_signature_path is None:
            payload_signature_path = str(
                _default_abis_resource_path(
                    technology=technology_normalized,
                    resource_name="sigmatrixRNAseq.txt" if technology_normalized == "rnaseq" else "sigmatrixMicro.txt",
                )
            )
        if technology_normalized == "microarray" and payload_target_path is None:
            payload_target_path = str(
                _default_abis_resource_path(
                    technology=technology_normalized,
                    resource_name="target.txt",
                )
            )

        result, _ = run_r_tool(
            tool="abis",
            payload={
                "mixture_path": str(mixture_path),
                "technology": technology_normalized,
                "signature_path": payload_signature_path,
                "target_path": "" if payload_target_path is None else payload_target_path,
                "maxit": 100,
            },
            rscript_command=rscript_command,
            output_path=output_path,
            log_path=log_path,
        )

    if rescale_percent:
        result = result / 100.0

    output_paths = {}
    if output_path is not None:
        output_paths["proportion"] = output_path
    if log_path is not None:
        output_paths["log"] = log_path

    return DeconvolutionResult(
        tool="abis",
        proportion=result,
        metadata={
            "transform": transform,
            "technology": technology_normalized,
            "signature_path": payload_signature_path,
            "target_path": "" if payload_target_path is None else payload_target_path,
            "native_output_scale": "percent",
            "output_scale": "fraction" if rescale_percent else "percent",
            "gene_match": "native_intersection",
            "r_env": "" if r_env is None else r_env,
        },
        output_paths=output_paths,
    )


def run_estimate(
    mixture: pd.DataFrame,
    transform: str = "cpm",
    platform: str = "illumina",
    gene_lengths: pd.Series | None = None,
    r_env: str | None = None,
    output_path: str | None = None,
    log_path: str | None = None,
) -> DeconvolutionResult:
    """Run ESTIMATE through the original R-Forge `estimate` package."""

    validate_expression_matrix(mixture, "mixture")
    mixture_norm = normalize_gene_names(mixture)
    mixture_transformed = apply_transform(mixture_norm, transform=transform, gene_lengths=gene_lengths)

    rscript_command = None
    if r_env is not None:
        rscript_command = ["mamba", "run", "-n", r_env, "Rscript"]

    with tempfile.TemporaryDirectory(prefix="foli_estimate_") as workdir:
        workdir_path = Path(workdir)
        mixture_path = workdir_path / "mixture.tsv"
        write_expression_matrix(mixture_transformed, str(mixture_path))
        result, _ = run_r_tool(
            tool="estimate",
            payload={
                "mixture_path": str(mixture_path),
                "platform": platform,
            },
            rscript_command=rscript_command,
            output_path=output_path,
            log_path=log_path,
        )

    output_paths = {}
    if output_path is not None:
        output_paths["score"] = output_path
    if log_path is not None:
        output_paths["log"] = log_path

    return DeconvolutionResult(
        tool="estimate",
        score=result,
        metadata={
            "transform": transform,
            "platform": platform,
            "r_env": "" if r_env is None else r_env,
            "output_semantics": "ESTIMATE stromal, immune, ESTIMATE, and optional tumor-purity scores; not cell-type proportions",
        },
        output_paths=output_paths,
    )


def run_consensus_tme(
    mixture: pd.DataFrame,
    transform: str = "cpm",
    cancer_type: str | None = None,
    stat_method: str = "ssgsea",
    sing_score_disp: bool = False,
    immune_score: bool = True,
    exclude_cells: list[str] | None = None,
    parallel_size: int = 0,
) -> DeconvolutionResult:
    """Run CONSENSUS-TME through its native R API."""

    validate_expression_matrix(mixture, "mixture")
    mixture_norm = normalize_gene_names(mixture)
    mixture_transformed = apply_transform(mixture_norm, transform=transform)

    with tempfile.TemporaryDirectory(prefix="foli_consensus_tme_") as workdir:
        workdir_path = Path(workdir)
        mixture_path = workdir_path / "mixture.tsv"
        write_expression_matrix(mixture_transformed, str(mixture_path))
        result, _ = run_r_tool(
            tool="consensus_tme",
            payload={
                "mixture_path": str(mixture_path),
                "cancer_type": "" if cancer_type is None else cancer_type,
                "stat_method": stat_method,
                "sing_score_disp": sing_score_disp,
                "immune_score": immune_score,
                "exclude_cells": exclude_cells,
                "parallel_size": parallel_size,
            },
        )

    return DeconvolutionResult(
        tool="consensus_tme",
        score=result,
        metadata={
            "transform": transform,
            "cancer_type": "" if cancer_type is None else cancer_type,
            "stat_method": stat_method,
            "sing_score_disp": sing_score_disp,
            "immune_score": immune_score,
            "exclude_cell_count": 0 if exclude_cells is None else len(exclude_cells),
            "parallel_size": parallel_size,
        },
    )

def run_cdseq(
    mixture: pd.DataFrame,
    signature: pd.DataFrame | None = None,
    scrna_counts: pd.DataFrame | None = None,
    cell_types: pd.Series | None = None,
    mixture_transform: str = "raw",
    signature_transform: str = "raw",
    gene_lengths: pd.Series | None = None,
    cell_type_number: int | list[int] | None = None,
    beta: float | list[float] | None = 0.5,
    alpha: float = 5.0,
    mcmc_iterations: int = 700,
    dilution_factor: float = 1.0,
    gene_subset_size: int | None = None,
    block_number: int = 1,
    cpu_number: int | None = 1,
    verbose: bool = False,
    print_progress_msg_to_file: int = 0,
    r_env: str | None = None,
    r_output_path: str | None = None,
    extra_output_dir: str | None = None,
) -> DeconvolutionResult:
    """Run CDSeq complete deconvolution with optional reference-based component labels."""

    validate_expression_matrix(mixture, "mixture")
    rscript_command = None
    if r_env is not None:
        rscript_command = ["mamba", "run", "-n", r_env, "Rscript"]

    mixture_norm = normalize_gene_names(mixture)
    mixture_transformed = apply_transform(
        mixture_norm,
        transform=mixture_transform,
        gene_lengths=gene_lengths,
    )

    reference_matrix = None
    if signature is not None:
        validate_expression_matrix(signature, "signature")
        signature_norm = normalize_gene_names(signature)
        signature_transformed = apply_transform(
            signature_norm,
            transform=signature_transform,
            gene_lengths=gene_lengths,
        )
        mixture_transformed, reference_matrix = align_mixture_and_signature(
            mixture_transformed,
            signature_transformed,
        )
    if scrna_counts is not None:
        validate_expression_matrix(scrna_counts, "scrna_counts")
        if cell_types is None:
            raise ValueError("cell_types are required when scrna_counts is provided")
        if len(cell_types) != scrna_counts.shape[1]:
            raise ValueError("cell_types length must match number of single-cell columns")
        scrna_norm = normalize_gene_names(scrna_counts)
        mixture_transformed, scrna_aligned = align_mixture_and_signature(
            mixture_transformed,
            scrna_norm,
        )
        aligned_cell_types = pd.Series(cell_types.to_numpy(), index=scrna_norm.columns).loc[scrna_aligned.columns]
        reference_matrix = scrna_aligned.T.groupby(aligned_cell_types).sum().T
        reference_matrix = apply_transform(
            reference_matrix,
            transform=signature_transform,
            gene_lengths=gene_lengths,
        )

    if cell_type_number is None:
        if reference_matrix is None:
            raise ValueError("cell_type_number is required when no signature or scRNA reference is provided")
        cell_type_number = int(reference_matrix.shape[1])

    with tempfile.TemporaryDirectory(prefix="foli_cdseq_") as workdir:
        workdir_path = Path(workdir)
        mixture_path = workdir_path / "mixture.tsv"
        write_expression_matrix(mixture_transformed, str(mixture_path))

        payload: dict[str, object] = {
            "mixture_path": str(mixture_path),
            "cell_type_number": [cell_type_number] if isinstance(cell_type_number, int) else cell_type_number,
            "beta": beta,
            "alpha": alpha,
            "mcmc_iterations": mcmc_iterations,
            "dilution_factor": dilution_factor,
            "gene_subset_size": gene_subset_size,
            "block_number": block_number,
            "cpu_number": cpu_number,
            "verbose": verbose,
            "print_progress_msg_to_file": print_progress_msg_to_file,
            "extra_output_dir": extra_output_dir,
        }
        if reference_matrix is not None:
            reference_path = workdir_path / "reference_gep.tsv"
            write_expression_matrix(reference_matrix, str(reference_path))
            payload["reference_gep_path"] = str(reference_path)
        if gene_lengths is not None:
            gene_lengths_path = workdir_path / "gene_lengths.tsv"
            gene_lengths.loc[mixture_transformed.index].to_frame(name="gene_length").to_csv(
                gene_lengths_path,
                sep="\t",
                index=False,
            )
            payload["gene_lengths_path"] = str(gene_lengths_path)

        result, _ = run_r_tool(
            tool="cdseq",
            payload=payload,
            rscript_command=rscript_command,
            output_path=r_output_path,
        )

    result_kwargs = {"component": result} if reference_matrix is None else {"proportion": result}
    return DeconvolutionResult(
        tool="cdseq",
        **result_kwargs,
        metadata={
            "mixture_transform": mixture_transform,
            "signature_transform": signature_transform,
            "reference_mode": "none" if reference_matrix is None else "signature_or_scrna",
            "cell_type_number": ",".join(str(value) for value in payload["cell_type_number"]),
            "mcmc_iterations": mcmc_iterations,
            "dilution_factor": dilution_factor,
            "gene_subset_size": 0 if gene_subset_size is None else gene_subset_size,
            "block_number": block_number,
            "cpu_number": 0 if cpu_number is None else cpu_number,
            "r_env": "" if r_env is None else r_env,
            "r_output_path": "" if r_output_path is None else r_output_path,
            "extra_output_dir": "" if extra_output_dir is None else extra_output_dir,
        },
    )


def run_xcell2(
    mixture: pd.DataFrame,
    xcell2_object_path: str | None = None,
    scrna_counts: pd.DataFrame | None = None,
    cell_types: pd.Series | None = None,
    dataset_labels: pd.Series | None = None,
    ontology_terms: pd.Series | None = None,
    xcell2_reference_output_path: str | None = None,
    transform: str = "cpm",
    scrna_transform: str = "raw",
    ref_type: str = "sc",
    use_ontology: bool = False,
    lineage_file: str | None = None,
    return_signatures: bool = False,
    training_use_spillover: bool = True,
    min_shared_genes: float = 0.2,
    raw_scores: bool = False,
    spillover: bool = True,
    spillover_alpha: float = 0.5,
    min_pb_cells: int = 30,
    min_pb_samples: int = 10,
    min_sc_genes: int = 10000,
    xcell2_workers: int = 1,
    seed: int = 427,
) -> DeconvolutionResult:
    """Run xCell2 against a trained object or train one from scRNA-seq first."""

    if xcell2_object_path is None:
        with tempfile.TemporaryDirectory(prefix="foli_xcell2_reference_") as workdir:
            workdir_path = Path(workdir)
            if xcell2_reference_output_path is None:
                object_path = workdir_path / "xcell2_reference.rds"
                reference_persisted = False
            else:
                object_path = Path(xcell2_reference_output_path)
                reference_persisted = True
            training_result = train_xcell2_reference(
                scrna_counts=scrna_counts,
                cell_types=cell_types,
                output_path=str(object_path),
                dataset_labels=dataset_labels,
                ontology_terms=ontology_terms,
                scrna_transform=scrna_transform,
                ref_type=ref_type,
                use_ontology=use_ontology,
                lineage_file=lineage_file,
                return_signatures=return_signatures,
                use_spillover=training_use_spillover,
                spillover_alpha=spillover_alpha,
                min_pb_cells=min_pb_cells,
                min_pb_samples=min_pb_samples,
                min_sc_genes=min_sc_genes,
                xcell2_workers=xcell2_workers,
                seed=seed,
            )
            result = run_xcell2(
                mixture=mixture,
                xcell2_object_path=str(object_path),
                transform=transform,
                min_shared_genes=min_shared_genes,
                raw_scores=True if return_signatures else raw_scores,
                spillover=False if return_signatures else spillover,
                spillover_alpha=spillover_alpha,
                xcell2_workers=xcell2_workers,
            )
            for key, value in training_result.metadata.items():
                result.metadata[f"training_{key}"] = value
            result.metadata["reference_persisted"] = reference_persisted
            return result

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
                "xcell2_workers": xcell2_workers,
            },
        )

    return DeconvolutionResult(
        tool="xcell2",
        score=result,
        metadata={
            "transform": transform,
            "min_shared_genes": min_shared_genes,
            "raw_scores": raw_scores,
            "spillover": spillover,
            "spillover_alpha": spillover_alpha,
            "xcell2_workers": xcell2_workers,
        },
    )


def train_xcell2_reference(
    scrna_counts: pd.DataFrame,
    cell_types: pd.Series,
    output_path: str,
    dataset_labels: pd.Series | None = None,
    ontology_terms: pd.Series | None = None,
    scrna_transform: str = "raw",
    ref_type: str = "sc",
    use_ontology: bool = False,
    lineage_file: str | None = None,
    return_signatures: bool = False,
    use_spillover: bool = True,
    spillover_alpha: float = 0.5,
    min_pb_cells: int = 30,
    min_pb_samples: int = 10,
    min_sc_genes: int = 10000,
    xcell2_workers: int = 1,
    seed: int = 427,
    log_path: str | None = None,
) -> ReferenceTrainingResult:
    """Train and save a custom xCell2 reference object."""

    validate_expression_matrix(scrna_counts, "scrna_counts")
    scrna_norm = normalize_gene_names(scrna_counts)
    scrna_transformed = apply_transform(scrna_norm, transform=scrna_transform)
    output = Path(output_path)
    output.parent.mkdir(parents=True, exist_ok=True)

    cell_types_aligned = cell_types.loc[scrna_transformed.columns].astype(str)
    if dataset_labels is None:
        datasets = pd.Series("reference", index=scrna_transformed.columns)
    else:
        datasets = dataset_labels.loc[scrna_transformed.columns].astype(str)
    if ontology_terms is None:
        ontologies = pd.Series(pd.NA, index=scrna_transformed.columns)
    else:
        ontologies = ontology_terms.loc[scrna_transformed.columns]

    labels = pd.DataFrame(
        {
            "sample": scrna_transformed.columns,
            "label": cell_types_aligned.to_numpy(),
            "dataset": datasets.to_numpy(),
            "ont": ontologies.to_numpy(),
        }
    )

    with tempfile.TemporaryDirectory(prefix="foli_xcell2_train_") as workdir:
        workdir_path = Path(workdir)
        scrna_path = workdir_path / "scrna.tsv"
        labels_path = workdir_path / "labels.tsv"
        write_expression_matrix(scrna_transformed, str(scrna_path))
        labels.to_csv(labels_path, sep="\t", index=False)
        run_r_tool(
            tool="xcell2_train",
            payload={
                "scrna_path": str(scrna_path),
                "labels_path": str(labels_path),
                "xcell2_object_path": str(output),
                "ref_type": ref_type,
                "use_ontology": use_ontology,
                "lineage_file": "" if lineage_file is None else lineage_file,
                "return_signatures": return_signatures,
                "use_spillover": use_spillover,
                "spillover_alpha": spillover_alpha,
                "min_pb_cells": min_pb_cells,
                "min_pb_samples": min_pb_samples,
                "min_sc_genes": min_sc_genes,
                "xcell2_workers": xcell2_workers,
                "seed": seed,
            },
            log_path=log_path,
        )

    return ReferenceTrainingResult(
        tool="xcell2",
        rds_path=str(output),
        metadata={
            "reference_path": str(output),
            "reference_mode": "signature_only" if return_signatures else "full",
            "raw_scores": bool(return_signatures),
            "spillover": bool((not return_signatures) and use_spillover),
            "return_signatures": bool(return_signatures),
            "scrna_transform": scrna_transform,
            "ref_type": ref_type,
            "use_ontology": bool(use_ontology),
            "n_cells": int(scrna_transformed.shape[1]),
            "n_genes": int(scrna_transformed.shape[0]),
            "n_cell_types": int(cell_types_aligned.nunique()),
            "min_pb_cells": int(min_pb_cells),
            "min_pb_samples": int(min_pb_samples),
            "min_sc_genes": int(min_sc_genes),
            "xcell2_workers": int(xcell2_workers),
            "seed": int(seed),
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
        proportion=result,
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
        proportion=result,
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
        proportion=result,
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
        proportion=result,
        metadata={
            "mixture_transform": mixture_transform,
            "signature_transform": signature_transform,
            "dwls_submethod": dwls_submethod,
        },
    )
