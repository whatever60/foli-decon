"""Wrappers for deconvolution tools that directly use scRNA-seq references."""

import json
import shutil
import tempfile
from pathlib import Path

import pandas as pd

from foli_decon.io import read_sample_by_celltype, write_expression_matrix, write_vector
from foli_decon.models import DeconvolutionResult
from foli_decon.preprocess import (
    align_mixture_and_scrna,
    apply_transform,
    intersect_gene_sets,
    normalize_gene_names,
    validate_expression_matrix,
)
from foli_decon.python_runner import run_command
from foli_decon.r_runner import run_r_tool


def _validate_scrna_metadata(
    scrna_counts: pd.DataFrame,
    cell_types: pd.Series,
    batch_ids: pd.Series,
) -> None:
    """Validate that scRNA metadata vectors align to cell columns."""

    if len(cell_types) != scrna_counts.shape[1]:
        raise ValueError("cell_types length must match number of scRNA columns")
    if len(batch_ids) != scrna_counts.shape[1]:
        raise ValueError("batch_ids length must match number of scRNA columns")


def run_music(
    mixture: pd.DataFrame,
    scrna_counts: pd.DataFrame | None = None,
    cell_types: pd.Series | None = None,
    batch_ids: pd.Series | None = None,
    mixture_transform: str = "raw",
    scrna_transform: str = "raw",
    select_ct: list[str] | None = None,
    scrna_sce_rds_path: str | None = None,
    cell_type_column: str = "cellType",
    batch_id_column: str = "sampleID",
    r_env: str | None = None,
    r_output_path: str | None = None,
    extra_output_dir: str | None = None,
) -> DeconvolutionResult:
    """Run MuSiC using bulk mixtures plus single-cell reference counts."""

    validate_expression_matrix(mixture, "mixture")
    rscript_command = None
    if r_env is not None:
        rscript_command = ["mamba", "run", "-n", r_env, "Rscript"]

    if scrna_sce_rds_path is not None:
        if scrna_transform != "raw":
            raise ValueError("scrna_transform must be 'raw' when scrna_sce_rds_path is used")
        mixture_norm = normalize_gene_names(mixture)
        mixture_transformed = apply_transform(mixture_norm, transform=mixture_transform)

        with tempfile.TemporaryDirectory(prefix="foli_music_sce_") as workdir:
            workdir_path = Path(workdir)
            mixture_path = workdir_path / "mixture.tsv"
            write_expression_matrix(mixture_transformed, str(mixture_path))
            result, _ = run_r_tool(
                tool="music_sce",
                payload={
                    "mixture_path": str(mixture_path),
                    "scrna_sce_rds_path": scrna_sce_rds_path,
                    "cell_type_column": cell_type_column,
                    "batch_id_column": batch_id_column,
                    "select_ct": select_ct,
                    "extra_output_dir": extra_output_dir,
                },
                rscript_command=rscript_command,
                output_path=r_output_path,
            )

        return DeconvolutionResult(
            tool="music",
            proportion=result,
            metadata={
                "mixture_transform": mixture_transform,
                "scrna_transform": scrna_transform,
                "select_ct_count": 0 if select_ct is None else len(select_ct),
                "reference_mode": "sce_rds",
                "cell_type_column": cell_type_column,
                "batch_id_column": batch_id_column,
                "r_env": "" if r_env is None else r_env,
                "r_output_path": "" if r_output_path is None else r_output_path,
                "extra_output_dir": "" if extra_output_dir is None else extra_output_dir,
            },
        )

    if scrna_counts is None or cell_types is None or batch_ids is None:
        raise ValueError("Provide scrna_counts, cell_types, and batch_ids, or provide scrna_sce_rds_path")

    validate_expression_matrix(scrna_counts, "scrna_counts")
    _validate_scrna_metadata(scrna_counts, cell_types, batch_ids)

    mixture_norm = normalize_gene_names(mixture)
    scrna_norm = normalize_gene_names(scrna_counts)
    mixture_aligned, scrna_aligned = align_mixture_and_scrna(mixture_norm, scrna_norm)

    mixture_transformed = apply_transform(mixture_aligned, transform=mixture_transform)
    scrna_transformed = apply_transform(scrna_aligned, transform=scrna_transform)

    with tempfile.TemporaryDirectory(prefix="foli_music_") as workdir:
        workdir_path = Path(workdir)
        mixture_path = workdir_path / "mixture.tsv"
        scrna_path = workdir_path / "scrna.tsv"
        cell_types_path = workdir_path / "cell_types.tsv"
        batch_ids_path = workdir_path / "batch_ids.tsv"

        write_expression_matrix(mixture_transformed, str(mixture_path))
        write_expression_matrix(scrna_transformed, str(scrna_path))
        write_vector(pd.Series(cell_types.to_numpy()), str(cell_types_path), "cell_type")
        write_vector(pd.Series(batch_ids.to_numpy()), str(batch_ids_path), "batch_id")

        result, _ = run_r_tool(
            tool="music",
            payload={
                "mixture_path": str(mixture_path),
                "scrna_path": str(scrna_path),
                "cell_types_path": str(cell_types_path),
                "batch_ids_path": str(batch_ids_path),
                "select_ct": select_ct,
            },
            rscript_command=rscript_command,
            output_path=r_output_path,
        )

    return DeconvolutionResult(
        tool="music",
        proportion=result,
        metadata={
            "mixture_transform": mixture_transform,
            "scrna_transform": scrna_transform,
            "select_ct_count": 0 if select_ct is None else len(select_ct),
            "reference_mode": "matrix",
            "r_env": "" if r_env is None else r_env,
            "r_output_path": "" if r_output_path is None else r_output_path,
        },
    )


def run_music2(
    control_mixture: pd.DataFrame,
    case_mixture: pd.DataFrame,
    scrna_counts: pd.DataFrame | None = None,
    cell_types: pd.Series | None = None,
    batch_ids: pd.Series | None = None,
    control_transform: str = "raw",
    case_transform: str = "raw",
    scrna_transform: str = "raw",
    select_ct: list[str] | None = None,
    scrna_sce_rds_path: str | None = None,
    cell_type_column: str = "cellType",
    batch_id_column: str = "sampleID",
    r_env: str | None = None,
    r_output_path: str | None = None,
    extra_output_dir: str | None = None,
    expr_low: float = 20.0,
    prop_r: float = 0.1,
    eps_c: float = 0.05,
    eps_r: float = 0.01,
    n_resample: int = 20,
    sample_prop: float = 0.5,
    cutoff_expr: float = 0.05,
    cutoff_fc: float = 2.0,
    cutoff_c: float = 0.05,
    cutoff_r: float = 0.01,
    maxiter: int = 200,
    ct_cov: bool = False,
    centered: bool = False,
    normalize: bool = False,
) -> DeconvolutionResult:
    """Run MuSiC2 using control/case bulk mixtures plus a single-cell reference."""

    validate_expression_matrix(control_mixture, "control_mixture")
    validate_expression_matrix(case_mixture, "case_mixture")
    rscript_command = None
    if r_env is not None:
        rscript_command = ["mamba", "run", "-n", r_env, "Rscript"]

    control_norm = normalize_gene_names(control_mixture)
    case_norm = normalize_gene_names(case_mixture)
    control_transformed = apply_transform(control_norm, transform=control_transform)
    case_transformed = apply_transform(case_norm, transform=case_transform)

    payload_base: dict[str, object] = {
        "select_ct": select_ct,
        "extra_output_dir": extra_output_dir,
        "expr_low": expr_low,
        "prop_r": prop_r,
        "eps_c": eps_c,
        "eps_r": eps_r,
        "n_resample": n_resample,
        "sample_prop": sample_prop,
        "cutoff_expr": cutoff_expr,
        "cutoff_fc": cutoff_fc,
        "cutoff_c": cutoff_c,
        "cutoff_r": cutoff_r,
        "maxiter": maxiter,
        "ct_cov": ct_cov,
        "centered": centered,
        "normalize": normalize,
    }

    if scrna_sce_rds_path is not None:
        if scrna_transform != "raw":
            raise ValueError("scrna_transform must be 'raw' when scrna_sce_rds_path is used")
        control_aligned, case_aligned = intersect_gene_sets([control_transformed, case_transformed])

        with tempfile.TemporaryDirectory(prefix="foli_music2_sce_") as workdir:
            workdir_path = Path(workdir)
            control_path = workdir_path / "control_mixture.tsv"
            case_path = workdir_path / "case_mixture.tsv"
            write_expression_matrix(control_aligned, str(control_path))
            write_expression_matrix(case_aligned, str(case_path))
            result, _ = run_r_tool(
                tool="music2_sce",
                payload={
                    **payload_base,
                    "control_mixture_path": str(control_path),
                    "case_mixture_path": str(case_path),
                    "scrna_sce_rds_path": scrna_sce_rds_path,
                    "cell_type_column": cell_type_column,
                    "batch_id_column": batch_id_column,
                },
                rscript_command=rscript_command,
                output_path=r_output_path,
            )

        return DeconvolutionResult(
            tool="music2",
            proportion=result,
            metadata={
                "control_transform": control_transform,
                "case_transform": case_transform,
                "scrna_transform": scrna_transform,
                "select_ct_count": 0 if select_ct is None else len(select_ct),
                "reference_mode": "sce_rds",
                "cell_type_column": cell_type_column,
                "batch_id_column": batch_id_column,
                "r_env": "" if r_env is None else r_env,
                "r_output_path": "" if r_output_path is None else r_output_path,
                "extra_output_dir": "" if extra_output_dir is None else extra_output_dir,
                "n_resample": n_resample,
                "sample_prop": sample_prop,
                "cutoff_fc": cutoff_fc,
                "maxiter": maxiter,
            },
        )

    if scrna_counts is None or cell_types is None or batch_ids is None:
        raise ValueError("Provide scrna_counts, cell_types, and batch_ids, or provide scrna_sce_rds_path")

    validate_expression_matrix(scrna_counts, "scrna_counts")
    _validate_scrna_metadata(scrna_counts, cell_types, batch_ids)
    scrna_norm = normalize_gene_names(scrna_counts)
    scrna_transformed = apply_transform(scrna_norm, transform=scrna_transform)
    control_aligned, case_aligned, scrna_aligned = intersect_gene_sets(
        [control_transformed, case_transformed, scrna_transformed]
    )

    with tempfile.TemporaryDirectory(prefix="foli_music2_") as workdir:
        workdir_path = Path(workdir)
        control_path = workdir_path / "control_mixture.tsv"
        case_path = workdir_path / "case_mixture.tsv"
        scrna_path = workdir_path / "scrna.tsv"
        cell_types_path = workdir_path / "cell_types.tsv"
        batch_ids_path = workdir_path / "batch_ids.tsv"

        write_expression_matrix(control_aligned, str(control_path))
        write_expression_matrix(case_aligned, str(case_path))
        write_expression_matrix(scrna_aligned, str(scrna_path))
        write_vector(pd.Series(cell_types.to_numpy()), str(cell_types_path), "cell_type")
        write_vector(pd.Series(batch_ids.to_numpy()), str(batch_ids_path), "batch_id")

        result, _ = run_r_tool(
            tool="music2",
            payload={
                **payload_base,
                "control_mixture_path": str(control_path),
                "case_mixture_path": str(case_path),
                "scrna_path": str(scrna_path),
                "cell_types_path": str(cell_types_path),
                "batch_ids_path": str(batch_ids_path),
            },
            rscript_command=rscript_command,
            output_path=r_output_path,
        )

    return DeconvolutionResult(
        tool="music2",
        proportion=result,
        metadata={
            "control_transform": control_transform,
            "case_transform": case_transform,
            "scrna_transform": scrna_transform,
            "select_ct_count": 0 if select_ct is None else len(select_ct),
            "reference_mode": "matrix",
            "r_env": "" if r_env is None else r_env,
            "r_output_path": "" if r_output_path is None else r_output_path,
            "extra_output_dir": "" if extra_output_dir is None else extra_output_dir,
            "n_resample": n_resample,
            "sample_prop": sample_prop,
            "cutoff_fc": cutoff_fc,
            "maxiter": maxiter,
        },
    )


def run_scdc(
    mixture: pd.DataFrame,
    scrna_counts: pd.DataFrame,
    cell_types: pd.Series,
    batch_ids: pd.Series,
    mixture_transform: str = "raw",
    scrna_transform: str = "raw",
    ct_sub: list[str] | None = None,
    iter_max: int = 1000,
    nu: float = 1e-4,
    epsilon: float = 0.001,
    weight_basis: bool = True,
    transform_bisque: bool = False,
    r_env: str | None = None,
    r_output_path: str | None = None,
    extra_output_dir: str | None = None,
) -> DeconvolutionResult:
    """Run SCDC with bulk mixtures plus single-cell reference counts."""

    validate_expression_matrix(mixture, "mixture")
    validate_expression_matrix(scrna_counts, "scrna_counts")
    _validate_scrna_metadata(scrna_counts, cell_types, batch_ids)
    rscript_command = None
    if r_env is not None:
        rscript_command = ["mamba", "run", "-n", r_env, "Rscript"]

    mixture_norm = normalize_gene_names(mixture)
    scrna_norm = normalize_gene_names(scrna_counts)
    mixture_transformed = apply_transform(mixture_norm, transform=mixture_transform)
    scrna_transformed = apply_transform(scrna_norm, transform=scrna_transform)
    mixture_aligned, scrna_aligned = intersect_gene_sets([mixture_transformed, scrna_transformed])

    with tempfile.TemporaryDirectory(prefix="foli_scdc_") as workdir:
        workdir_path = Path(workdir)
        mixture_path = workdir_path / "mixture.tsv"
        scrna_path = workdir_path / "scrna.tsv"
        cell_types_path = workdir_path / "cell_types.tsv"
        batch_ids_path = workdir_path / "batch_ids.tsv"

        write_expression_matrix(mixture_aligned, str(mixture_path))
        write_expression_matrix(scrna_aligned, str(scrna_path))
        write_vector(pd.Series(cell_types.to_numpy()), str(cell_types_path), "cell_type")
        write_vector(pd.Series(batch_ids.to_numpy()), str(batch_ids_path), "batch_id")

        result, _ = run_r_tool(
            tool="scdc",
            payload={
                "mixture_path": str(mixture_path),
                "scrna_path": str(scrna_path),
                "cell_types_path": str(cell_types_path),
                "batch_ids_path": str(batch_ids_path),
                "ct_sub": ct_sub,
                "iter_max": iter_max,
                "nu": nu,
                "epsilon": epsilon,
                "weight_basis": weight_basis,
                "transform_bisque": transform_bisque,
                "extra_output_dir": extra_output_dir,
            },
            rscript_command=rscript_command,
            output_path=r_output_path,
        )

    return DeconvolutionResult(
        tool="scdc",
        proportion=result,
        metadata={
            "mixture_transform": mixture_transform,
            "scrna_transform": scrna_transform,
            "ct_sub_count": 0 if ct_sub is None else len(ct_sub),
            "iter_max": iter_max,
            "nu": nu,
            "epsilon": epsilon,
            "weight_basis": weight_basis,
            "transform_bisque": transform_bisque,
            "r_env": "" if r_env is None else r_env,
            "r_output_path": "" if r_output_path is None else r_output_path,
            "extra_output_dir": "" if extra_output_dir is None else extra_output_dir,
        },
    )


def run_mead(
    mixture: pd.DataFrame,
    scrna_counts: pd.DataFrame,
    cell_types: pd.Series,
    batch_ids: pd.Series,
    mixture_transform: str = "raw",
    scrna_transform: str = "raw",
    select_ct: list[str] | None = None,
    marker_genes: list[str] | None = None,
    gene_thresh: float = 0.0,
    max_count_quantile_celltype: float | None = None,
    max_count_quantile_indi: float | None = None,
    filter_gene: bool = False,
    hc_type: str = "hc3",
    centering_xy: bool = False,
    nfold: int = 10,
    groups: pd.Series | None = None,
    calc_var: bool = True,
    use_qp: bool = False,
    r_env: str | None = None,
    r_output_path: str | None = None,
    extra_output_dir: str | None = None,
) -> DeconvolutionResult:
    """Run MEAD with bulk mixtures plus a SingleCellExperiment-style scRNA reference."""

    validate_expression_matrix(mixture, "mixture")
    validate_expression_matrix(scrna_counts, "scrna_counts")
    _validate_scrna_metadata(scrna_counts, cell_types, batch_ids)
    if groups is not None and len(groups) != mixture.shape[1]:
        raise ValueError("groups length must match number of mixture samples")

    rscript_command = None
    if r_env is not None:
        rscript_command = ["mamba", "run", "-n", r_env, "Rscript"]

    mixture_norm = normalize_gene_names(mixture)
    scrna_norm = normalize_gene_names(scrna_counts)
    mixture_transformed = apply_transform(mixture_norm, transform=mixture_transform)
    scrna_transformed = apply_transform(scrna_norm, transform=scrna_transform)
    mixture_aligned, scrna_aligned = intersect_gene_sets([mixture_transformed, scrna_transformed])

    uncertainty = None
    with tempfile.TemporaryDirectory(prefix="foli_mead_") as workdir:
        workdir_path = Path(workdir)
        mixture_path = workdir_path / "mixture.tsv"
        scrna_path = workdir_path / "scrna.tsv"
        cell_types_path = workdir_path / "cell_types.tsv"
        batch_ids_path = workdir_path / "batch_ids.tsv"
        uncertainty_path = workdir_path / "mead_p_hat_se.tsv"

        write_expression_matrix(mixture_aligned, str(mixture_path))
        write_expression_matrix(scrna_aligned, str(scrna_path))
        write_vector(pd.Series(cell_types.to_numpy()), str(cell_types_path), "cell_type")
        write_vector(pd.Series(batch_ids.to_numpy()), str(batch_ids_path), "batch_id")

        payload: dict[str, object] = {
            "mixture_path": str(mixture_path),
            "scrna_path": str(scrna_path),
            "cell_types_path": str(cell_types_path),
            "batch_ids_path": str(batch_ids_path),
            "select_ct": select_ct,
            "marker_genes": marker_genes,
            "gene_thresh": gene_thresh,
            "max_count_quantile_celltype": max_count_quantile_celltype,
            "max_count_quantile_indi": max_count_quantile_indi,
            "filter_gene": filter_gene,
            "hc_type": hc_type,
            "centering_xy": centering_xy,
            "nfold": nfold,
            "calc_var": calc_var,
            "use_qp": use_qp,
            "uncertainty_path": str(uncertainty_path),
            "extra_output_dir": extra_output_dir,
        }
        if groups is not None:
            groups_path = workdir_path / "groups.tsv"
            write_vector(pd.Series(groups.to_numpy()), str(groups_path), "group")
            payload["groups_path"] = str(groups_path)

        result, _ = run_r_tool(
            tool="mead",
            payload=payload,
            rscript_command=rscript_command,
            output_path=r_output_path,
        )
        if calc_var:
            uncertainty = read_sample_by_celltype(str(uncertainty_path))

    output_paths = {}
    if r_output_path is not None:
        output_paths["proportion"] = r_output_path
    if extra_output_dir is not None:
        output_paths["extra_output_dir"] = extra_output_dir

    return DeconvolutionResult(
        tool="mead",
        proportion=result,
        uncertainty=uncertainty,
        metadata={
            "mixture_transform": mixture_transform,
            "scrna_transform": scrna_transform,
            "select_ct_count": 0 if select_ct is None else len(select_ct),
            "marker_gene_count": 0 if marker_genes is None else len(marker_genes),
            "gene_thresh": gene_thresh,
            "filter_gene": filter_gene,
            "hc_type": hc_type,
            "centering_xy": centering_xy,
            "nfold": nfold,
            "groups_provided": groups is not None,
            "calc_var": calc_var,
            "use_qp": use_qp,
            "r_env": "" if r_env is None else r_env,
            "r_output_path": "" if r_output_path is None else r_output_path,
            "extra_output_dir": "" if extra_output_dir is None else extra_output_dir,
            "output_semantics": "MEAD p_hat cell-type proportions; uncertainty stores p_hat standard errors",
        },
        output_paths=output_paths,
    )


def run_instaprism(
    mixture: pd.DataFrame,
    scrna_counts: pd.DataFrame,
    cell_types: pd.Series,
    mixture_transform: str = "raw",
    scrna_transform: str = "raw",
    cell_states: pd.Series | None = None,
    r_env: str | None = None,
    r_output_path: str | None = None,
    extra_output_dir: str | None = None,
    write_z: bool = False,
) -> DeconvolutionResult:
    """Run InstaPrism with a reference prepared from single-cell counts."""

    validate_expression_matrix(mixture, "mixture")
    validate_expression_matrix(scrna_counts, "scrna_counts")
    if len(cell_types) != scrna_counts.shape[1]:
        raise ValueError("cell_types length must match number of scRNA columns")
    if cell_states is not None and len(cell_states) != scrna_counts.shape[1]:
        raise ValueError("cell_states length must match number of scRNA columns")

    rscript_command = None
    if r_env is not None:
        rscript_command = ["mamba", "run", "-n", r_env, "Rscript"]

    mixture_norm = normalize_gene_names(mixture)
    scrna_norm = normalize_gene_names(scrna_counts)
    mixture_transformed = apply_transform(mixture_norm, transform=mixture_transform)
    scrna_transformed = apply_transform(scrna_norm, transform=scrna_transform)
    mixture_aligned, scrna_aligned = intersect_gene_sets([mixture_transformed, scrna_transformed])

    with tempfile.TemporaryDirectory(prefix="foli_instaprism_") as workdir:
        workdir_path = Path(workdir)
        mixture_path = workdir_path / "mixture.tsv"
        scrna_path = workdir_path / "scrna.tsv"
        cell_types_path = workdir_path / "cell_types.tsv"

        write_expression_matrix(mixture_aligned, str(mixture_path))
        write_expression_matrix(scrna_aligned, str(scrna_path))
        write_vector(pd.Series(cell_types.to_numpy()), str(cell_types_path), "cell_type")

        payload: dict[str, object] = {
            "mixture_path": str(mixture_path),
            "scrna_path": str(scrna_path),
            "cell_types_path": str(cell_types_path),
            "extra_output_dir": extra_output_dir,
            "write_z": write_z,
        }
        if cell_states is not None:
            cell_states_path = workdir_path / "cell_states.tsv"
            write_vector(pd.Series(cell_states.to_numpy()), str(cell_states_path), "cell_state")
            payload["cell_states_path"] = str(cell_states_path)

        result, _ = run_r_tool(
            tool="instaprism",
            payload=payload,
            rscript_command=rscript_command,
            output_path=r_output_path,
        )

    return DeconvolutionResult(
        tool="instaprism",
        proportion=result,
        metadata={
            "mixture_transform": mixture_transform,
            "scrna_transform": scrna_transform,
            "cell_state_mode": "cell_types" if cell_states is None else "provided",
            "r_env": "" if r_env is None else r_env,
            "r_output_path": "" if r_output_path is None else r_output_path,
            "extra_output_dir": "" if extra_output_dir is None else extra_output_dir,
            "write_z": write_z,
        },
    )


def run_tape(
    mixture: pd.DataFrame,
    scrna_counts: pd.DataFrame,
    cell_types: pd.Series,
    mixture_transform: str = "raw",
    scrna_transform: str = "raw",
    python_env: str | None = None,
    variance_threshold: float = 0.98,
    scaler: str = "mms",
    datatype: str = "counts",
    gene_lengths_table_path: str | None = None,
    mode: str = "overall",
    adaptive: bool = True,
    sparse: bool = True,
    batch_size: int = 128,
    epochs: int = 128,
    seed: int = 0,
    output_path: str | None = None,
    signature_output_path: str | None = None,
    log_path: str | None = None,
) -> DeconvolutionResult:
    """Run TAPE with a single-cell reference through a CPU PyTorch sidecar."""

    validate_expression_matrix(mixture, "mixture")
    validate_expression_matrix(scrna_counts, "scrna_counts")
    if len(cell_types) != scrna_counts.shape[1]:
        raise ValueError("cell_types length must match number of scRNA columns")

    command_base = ["python"] if python_env is None else ["mamba", "run", "-n", python_env, "python"]
    runner_path = Path(__file__).resolve().parents[1] / "deep_scripts" / "run_tape.py"

    mixture_norm = normalize_gene_names(mixture)
    scrna_norm = normalize_gene_names(scrna_counts)
    mixture_transformed = apply_transform(mixture_norm, transform=mixture_transform)
    scrna_transformed = apply_transform(scrna_norm, transform=scrna_transform)
    mixture_aligned, scrna_aligned = intersect_gene_sets([mixture_transformed, scrna_transformed])
    aligned_cell_types = pd.Series(cell_types.to_numpy(), index=scrna_norm.columns).loc[scrna_aligned.columns]

    with tempfile.TemporaryDirectory(prefix="foli_tape_") as workdir:
        workdir_path = Path(workdir)
        mixture_path = workdir_path / "mixture.tsv"
        scrna_path = workdir_path / "scrna_reference.tsv"
        payload_path = workdir_path / "payload.json"
        resolved_output_path = (
            workdir_path / "tape_proportions.tsv" if output_path is None else Path(output_path).resolve()
        )
        resolved_signature_output_path = (
            "" if signature_output_path is None else str(Path(signature_output_path).resolve())
        )

        write_expression_matrix(mixture_aligned.T, str(mixture_path))
        scrna_for_tape = scrna_aligned.T.copy()
        scrna_for_tape.index = aligned_cell_types.astype(str).to_numpy()
        write_expression_matrix(scrna_for_tape, str(scrna_path))
        payload = {
            "mixture_path": str(mixture_path),
            "scrna_path": str(scrna_path),
            "output_path": str(resolved_output_path),
            "signature_output_path": resolved_signature_output_path,
            "variance_threshold": variance_threshold,
            "scaler": scaler,
            "datatype": datatype,
            "gene_lengths_table_path": "" if gene_lengths_table_path is None else gene_lengths_table_path,
            "mode": mode,
            "adaptive": adaptive,
            "sparse": sparse,
            "batch_size": batch_size,
            "epochs": epochs,
            "seed": seed,
        }
        payload_path.write_text(json.dumps(payload), encoding="utf-8")
        run_command(command_base + [str(runner_path), str(payload_path)], cwd=str(workdir_path), log_path=log_path)
        result = read_sample_by_celltype(str(resolved_output_path))

    output_paths = {}
    if output_path is not None:
        output_paths["proportion"] = output_path
    if signature_output_path is not None:
        output_paths["signature"] = signature_output_path

    return DeconvolutionResult(
        tool="tape",
        proportion=result,
        metadata={
            "mixture_transform": mixture_transform,
            "scrna_transform": scrna_transform,
            "python_env": "" if python_env is None else python_env,
            "variance_threshold": variance_threshold,
            "scaler": scaler,
            "datatype": datatype,
            "mode": mode,
            "adaptive": adaptive,
            "sparse": sparse,
            "batch_size": batch_size,
            "epochs": epochs,
            "seed": seed,
            "log_path": "" if log_path is None else log_path,
        },
        output_paths=output_paths,
    )


def run_scaden(
    mixture: pd.DataFrame,
    scrna_counts: pd.DataFrame | None = None,
    cell_types: pd.Series | None = None,
    mixture_transform: str = "raw",
    scrna_transform: str = "raw",
    python_env: str | None = None,
    model_dir: str | None = None,
    retrain: bool = True,
    n_training_samples: int = 1000,
    cells_per_sample: int = 100,
    train_steps: int = 5000,
    batch_size: int = 128,
    learning_rate: float = 0.0001,
    var_cutoff: float = 0.0,
    seed: int = 0,
    unknown_celltypes: list[str] | None = None,
    output_path: str | None = None,
    log_path: str | None = None,
) -> DeconvolutionResult:
    """Run Scaden through its simulate/process/train/predict command-line workflow."""

    validate_expression_matrix(mixture, "mixture")
    command_base = ["scaden"] if python_env is None else ["mamba", "run", "-n", python_env, "scaden"]

    mixture_norm = normalize_gene_names(mixture)
    mixture_transformed = apply_transform(mixture_norm, transform=mixture_transform)

    with tempfile.TemporaryDirectory(prefix="foli_scaden_") as workdir:
        workdir_path = Path(workdir)
        mixture_path = workdir_path / "mixture.tsv"
        prediction_path = workdir_path / "scaden_predictions.txt" if output_path is None else Path(output_path)

        if retrain:
            if scrna_counts is None or cell_types is None:
                raise ValueError("scrna_counts and cell_types are required when retrain=True")
            validate_expression_matrix(scrna_counts, "scrna_counts")
            if len(cell_types) != scrna_counts.shape[1]:
                raise ValueError("cell_types length must match number of scRNA columns")

            scrna_norm = normalize_gene_names(scrna_counts)
            scrna_transformed = apply_transform(scrna_norm, transform=scrna_transform)
            mixture_transformed, scrna_aligned = intersect_gene_sets([mixture_transformed, scrna_transformed])
            aligned_cell_types = pd.Series(cell_types.to_numpy(), index=scrna_norm.columns).loc[scrna_aligned.columns]

            scaden_reference = scrna_aligned.T.copy()
            scaden_reference.index = range(scaden_reference.shape[0])
            scaden_reference.to_csv(workdir_path / "reference_counts.txt", sep="\t")
            aligned_cell_types.to_frame(name="Celltype").to_csv(
                workdir_path / "reference_celltypes.txt",
                sep="\t",
                index=False,
            )

            simulation_dir = workdir_path / "simulation"
            simulation_dir.mkdir()
            unknown_values = ["__foli_no_unknown_celltype__"] if unknown_celltypes is None else unknown_celltypes
            simulate_command = command_base + [
                "simulate",
                "--out",
                str(simulation_dir),
                "--data",
                str(workdir_path),
                "--cells",
                str(cells_per_sample),
                "--n_samples",
                str(n_training_samples),
                "--pattern",
                "*_counts.txt",
                "--prefix",
                str(simulation_dir / "training"),
                "--data-format",
                "txt",
            ]
            for unknown_celltype in unknown_values:
                simulate_command.extend(["--unknown", unknown_celltype])

            write_expression_matrix(mixture_transformed, str(mixture_path))
            processed_path = workdir_path / "processed.h5ad"
            resolved_model_dir = workdir_path / "model" if model_dir is None else Path(model_dir)
            resolved_model_dir.parent.mkdir(parents=True, exist_ok=True)

            run_command(simulate_command, log_path=log_path)
            run_command(
                command_base
                + [
                    "process",
                    str(simulation_dir / "training.h5ad"),
                    str(mixture_path),
                    "--processed_path",
                    str(processed_path),
                    "--var_cutoff",
                    str(var_cutoff),
                ],
                log_path=log_path,
            )
            run_command(
                command_base
                + [
                    "train",
                    str(processed_path),
                    "--model_dir",
                    str(resolved_model_dir),
                    "--batch_size",
                    str(batch_size),
                    "--learning_rate",
                    str(learning_rate),
                    "--steps",
                    str(train_steps),
                    "--seed",
                    str(seed),
                ],
                log_path=log_path,
            )
        else:
            if model_dir is None:
                raise ValueError("model_dir is required when retrain=False")
            resolved_model_dir = Path(model_dir)
            write_expression_matrix(mixture_transformed, str(mixture_path))

        prediction_path.parent.mkdir(parents=True, exist_ok=True)
        run_command(
            command_base
            + [
                "predict",
                str(mixture_path),
                "--model_dir",
                str(resolved_model_dir),
                "--outname",
                str(prediction_path),
                "--seed",
                str(seed),
            ],
            log_path=log_path,
        )
        result = read_sample_by_celltype(str(prediction_path))

    return DeconvolutionResult(
        tool="scaden",
        proportion=result,
        metadata={
            "mixture_transform": mixture_transform,
            "scrna_transform": scrna_transform,
            "python_env": "" if python_env is None else python_env,
            "model_dir": "" if model_dir is None else model_dir,
            "retrain": retrain,
            "n_training_samples": n_training_samples,
            "cells_per_sample": cells_per_sample,
            "train_steps": train_steps,
            "batch_size": batch_size,
            "learning_rate": learning_rate,
            "var_cutoff": var_cutoff,
            "seed": seed,
            "output_path": "" if output_path is None else output_path,
            "log_path": "" if log_path is None else log_path,
        },
    )


def run_dissect(
    mixture: pd.DataFrame,
    scrna_counts: pd.DataFrame,
    cell_types: pd.Series,
    batch_ids: pd.Series | None = None,
    mixture_transform: str = "cpm",
    scrna_transform: str = "raw",
    python_env: str | None = None,
    n_training_samples: int = 1000,
    cells_per_sample: int = 500,
    prop_sparse: float = 0.5,
    train_steps: int = 5000,
    batch_size: int = 64,
    learning_rate: float = 1e-5,
    n_models: int = 5,
    var_cutoff: float = 0.0,
    alpha_range: list[float] | None = None,
    test_in_mix: int | None = None,
    normalize_simulated: str | None = "cpm",
    normalization_per_batch: str | None = "log1p-MinMax",
    mix: str = "srm",
    min_genes: int = 0,
    min_cells: int = 0,
    mt_cutoff: float = 100.0,
    min_expr: float = 0.0,
    output_path: str | None = None,
    extra_output_dir: str | None = None,
    log_path: str | None = None,
) -> DeconvolutionResult:
    """Run DISSECT fraction estimation through a CPU TensorFlow sidecar environment."""

    validate_expression_matrix(mixture, "mixture")
    validate_expression_matrix(scrna_counts, "scrna_counts")
    if len(cell_types) != scrna_counts.shape[1]:
        raise ValueError("cell_types length must match number of scRNA columns")
    if batch_ids is not None and len(batch_ids) != scrna_counts.shape[1]:
        raise ValueError("batch_ids length must match number of scRNA columns")

    command_base = ["python"] if python_env is None else ["mamba", "run", "-n", python_env, "python"]
    runner_path = Path(__file__).resolve().parents[1] / "deep_scripts" / "run_dissect_cpu.py"

    mixture_norm = normalize_gene_names(mixture)
    scrna_norm = normalize_gene_names(scrna_counts)
    mixture_transformed = apply_transform(mixture_norm, transform=mixture_transform)
    scrna_transformed = apply_transform(scrna_norm, transform=scrna_transform)
    mixture_aligned, scrna_aligned = intersect_gene_sets([mixture_transformed, scrna_transformed])
    aligned_cell_types = pd.Series(cell_types.to_numpy(), index=scrna_norm.columns).loc[scrna_aligned.columns]

    previous_string_storage = pd.get_option("mode.string_storage")
    restore_string_storage = previous_string_storage != "python"
    if restore_string_storage:
        pd.set_option("mode.string_storage", "python")

    with tempfile.TemporaryDirectory(prefix="foli_dissect_") as workdir:
        workdir_path = Path(workdir)
        mixture_path = workdir_path / "mixture.tsv"
        scrna_matrix_path = workdir_path / "scrna_reference.tsv"
        scrna_metadata_path = workdir_path / "scrna_metadata.tsv"
        scrna_h5ad_path = workdir_path / "scrna_reference.h5ad"
        experiment_folder = workdir_path / "dissect_experiment"
        payload_path = workdir_path / "payload.json"

        write_expression_matrix(mixture_aligned, str(mixture_path))
        cell_ids = [str(value) for value in scrna_aligned.columns]
        metadata = pd.DataFrame(
            {
                "cell_id": cell_ids,
                "cell_type": [str(value) for value in aligned_cell_types.to_numpy()],
            }
        )
        batch_col = None
        if batch_ids is not None:
            aligned_batches = pd.Series(
                [str(v) for v in batch_ids.to_numpy()],
                index=[str(v) for v in batch_ids.index],
            ).loc[cell_ids]
            metadata["batch"] = [str(v) for v in aligned_batches.to_numpy()]
            batch_col = "batch"
        write_expression_matrix(scrna_aligned, str(scrna_matrix_path))
        metadata.to_csv(scrna_metadata_path, sep="\t", index=False)

        payload = {
            "mixture_path": str(mixture_path),
            "scrna_matrix_path": str(scrna_matrix_path),
            "scrna_metadata_path": str(scrna_metadata_path),
            "scrna_h5ad_path": str(scrna_h5ad_path),
            "experiment_folder": str(experiment_folder),
            "batch_col": batch_col,
            "n_training_samples": n_training_samples,
            "cells_per_sample": cells_per_sample,
            "prop_sparse": prop_sparse,
            "train_steps": train_steps,
            "batch_size": batch_size,
            "learning_rate": learning_rate,
            "n_models": n_models,
            "var_cutoff": var_cutoff,
            "alpha_range": [0.1, 0.9] if alpha_range is None else alpha_range,
            "test_in_mix": test_in_mix,
            "normalize_test": None,
            "normalize_simulated": normalize_simulated,
            "normalization_per_batch": normalization_per_batch,
            "mix": mix,
            "min_genes": min_genes,
            "min_cells": min_cells,
            "mt_cutoff": mt_cutoff,
            "min_expr": min_expr,
        }
        payload_path.write_text(json.dumps(payload), encoding="utf-8")
        run_command(command_base + [str(runner_path), str(payload_path)], log_path=log_path)

        native_result_path = experiment_folder / "dissect_fractions.txt"
        native_scores_path = experiment_folder / "dissect_scores.txt"
        result = read_sample_by_celltype(str(native_result_path))
        if output_path is not None:
            result.to_csv(output_path, sep="\t")
        if extra_output_dir is not None:
            extra_path = Path(extra_output_dir)
            extra_path.mkdir(parents=True, exist_ok=True)
            shutil.copyfile(native_result_path, extra_path / "dissect_fractions.txt")
            shutil.copyfile(native_scores_path, extra_path / "dissect_scores.txt")
            shutil.copyfile(experiment_folder / "main_config.json", extra_path / "main_config.json")

    if restore_string_storage:
        pd.set_option("mode.string_storage", previous_string_storage)

    return DeconvolutionResult(
        tool="dissect",
        proportion=result,
        metadata={
            "mixture_transform": mixture_transform,
            "scrna_transform": scrna_transform,
            "python_env": "" if python_env is None else python_env,
            "n_training_samples": n_training_samples,
            "cells_per_sample": cells_per_sample,
            "prop_sparse": prop_sparse,
            "train_steps": train_steps,
            "batch_size": batch_size,
            "learning_rate": learning_rate,
            "n_models": n_models,
            "var_cutoff": var_cutoff,
            "normalize_simulated": "" if normalize_simulated is None else normalize_simulated,
            "normalization_per_batch": "" if normalization_per_batch is None else normalization_per_batch,
            "mix": mix,
            "output_path": "" if output_path is None else output_path,
            "extra_output_dir": "" if extra_output_dir is None else extra_output_dir,
            "log_path": "" if log_path is None else log_path,
        },
    )


def run_blue(
    mixture: pd.DataFrame,
    scrna_counts: pd.DataFrame,
    cell_types: pd.Series,
    batch_ids: pd.Series | None = None,
    mixture_transform: str = "raw",
    scrna_transform: str = "raw",
    python_env: str | None = None,
    blue_repo_path: str | None = None,
    celltype_mapping: dict[str, list[str]] | None = None,
    input_genes: list[str] | None = None,
    output_genes: list[str] | None = None,
    samplenum_per_ct: int = 30,
    val_samplenum_per_patient: int = 10,
    n_training_samples: int = 3000,
    n_validation_samples: int = 300,
    cells_per_sample: int = 500,
    uniform_alpha: float = 1.0,
    dominant_alphas: list[float] | None = None,
    deg_fdr_rate: float = 0.1,
    deg_num_per_ct: int = 200,
    epochs: int = 50,
    batch_size: int = 64,
    learning_rate: float = 1e-5,
    weight_decay: float = 0.0,
    coeff_ct_gep: float = 1.0,
    coeff_prop: float = 100.0,
    sched_step: int = 20,
    sched_gamma: float = 0.5,
    preprocess_mode: str = "deg_first",
    split_mode: str = "pseudo",
    split_seed: int = 42,
    seed: int = 42,
    output_path: str | None = None,
    ctgep_output_path: str | None = None,
    extra_output_dir: str | None = None,
    log_path: str | None = None,
) -> DeconvolutionResult:
    """Run BLUE by staging inputs and executing the native BLUE pipeline on CPU."""

    validate_expression_matrix(mixture, "mixture")
    validate_expression_matrix(scrna_counts, "scrna_counts")
    if len(cell_types) != scrna_counts.shape[1]:
        raise ValueError("cell_types length must match number of scRNA columns")
    if batch_ids is not None and len(batch_ids) != scrna_counts.shape[1]:
        raise ValueError("batch_ids length must match number of scRNA columns")

    command_base = ["python"] if python_env is None else ["mamba", "run", "-n", python_env, "python"]
    runner_path = Path(__file__).resolve().parents[1] / "deep_scripts" / "run_blue_cpu.py"
    repo_root = Path(__file__).resolve().parents[3]
    resolved_blue_repo = repo_root / "resources" / "blue" / "BLUE" if blue_repo_path is None else Path(blue_repo_path)

    mixture_norm = normalize_gene_names(mixture)
    scrna_norm = normalize_gene_names(scrna_counts)
    mixture_transformed = apply_transform(mixture_norm, transform=mixture_transform)
    scrna_transformed = apply_transform(scrna_norm, transform=scrna_transform)
    mixture_aligned, scrna_aligned = intersect_gene_sets([mixture_transformed, scrna_transformed])
    aligned_cell_types = pd.Series(cell_types.to_numpy(), index=scrna_norm.columns).loc[scrna_aligned.columns].astype(str)
    if batch_ids is None:
        aligned_batch_ids = pd.Series(["reference"] * scrna_aligned.shape[1], index=scrna_aligned.columns)
    else:
        aligned_batch_ids = pd.Series(batch_ids.to_numpy(), index=scrna_norm.columns).loc[scrna_aligned.columns].astype(str)

    if celltype_mapping is None:
        resolved_mapping = {cell_type: [cell_type] for cell_type in sorted(aligned_cell_types.unique().tolist())}
    else:
        resolved_mapping = {
            str(coarse_type): [str(fine_type) for fine_type in fine_types]
            for coarse_type, fine_types in celltype_mapping.items()
        }
        mapped_fine_types: set[str] = set()
        for fine_types in resolved_mapping.values():
            mapped_fine_types.update(fine_types)
        missing_fine_types = sorted(set(aligned_cell_types.unique().tolist()) - mapped_fine_types)
        if missing_fine_types:
            raise ValueError(f"celltype_mapping does not cover these cell types: {missing_fine_types}")

    aligned_gene_set = set(mixture_aligned.index.astype(str))
    if input_genes is None:
        resolved_input_genes = [str(gene) for gene in mixture_aligned.index]
    else:
        resolved_input_genes = [str(gene).replace(".", "-") for gene in input_genes]
        missing_input_genes = sorted(set(resolved_input_genes) - aligned_gene_set)
        if missing_input_genes:
            raise ValueError(f"input_genes are not present in the shared mixture/scRNA genes: {missing_input_genes}")
    if output_genes is None:
        resolved_output_genes = resolved_input_genes
    else:
        resolved_output_genes = [str(gene).replace(".", "-") for gene in output_genes]
        missing_output_genes = sorted(set(resolved_output_genes) - aligned_gene_set)
        if missing_output_genes:
            raise ValueError(f"output_genes are not present in the shared mixture/scRNA genes: {missing_output_genes}")

    with tempfile.TemporaryDirectory(prefix="foli_blue_") as workdir:
        workdir_path = Path(workdir)
        mixture_path = workdir_path / "bulk.tsv"
        scrna_matrix_path = workdir_path / "scrna_reference.tsv"
        scrna_metadata_path = workdir_path / "scrna_metadata.tsv"
        input_gene_list_path = workdir_path / "blue_input_genes.txt"
        output_gene_list_path = workdir_path / "blue_output_genes.txt"
        payload_path = workdir_path / "payload.json"
        if extra_output_dir is None:
            resolved_output_path = workdir_path / "blue_proportions.tsv" if output_path is None else Path(output_path)
            resolved_ctgep_output_path = (
                workdir_path / "blue_predicted_ctGEP.h5ad" if ctgep_output_path is None else Path(ctgep_output_path)
            )
        else:
            extra_path = Path(extra_output_dir)
            extra_path.mkdir(parents=True, exist_ok=True)
            resolved_output_path = extra_path / "blue_proportions.tsv" if output_path is None else Path(output_path)
            resolved_ctgep_output_path = (
                extra_path / "blue_predicted_ctGEP.h5ad" if ctgep_output_path is None else Path(ctgep_output_path)
            )

        mixture_aligned.to_csv(mixture_path, sep="\t", index_label="gene_id")
        write_expression_matrix(scrna_aligned, str(scrna_matrix_path))
        pd.DataFrame(
            {
                "cell_id": [str(value) for value in scrna_aligned.columns],
                "cell_type": [str(value) for value in aligned_cell_types.to_numpy()],
                "library_id": [str(value) for value in aligned_batch_ids.to_numpy()],
            }
        ).to_csv(scrna_metadata_path, sep="\t", index=False)
        input_gene_list_path.write_text("\n".join(resolved_input_genes) + "\n", encoding="utf-8")
        output_gene_list_path.write_text("\n".join(resolved_output_genes) + "\n", encoding="utf-8")

        payload = {
            "blue_repo_path": str(resolved_blue_repo),
            "workdir": str(workdir_path / "blue_workdir"),
            "mixture_path": str(mixture_path),
            "scrna_matrix_path": str(scrna_matrix_path),
            "scrna_metadata_path": str(scrna_metadata_path),
            "input_gene_list_path": str(input_gene_list_path),
            "output_gene_list_path": str(output_gene_list_path),
            "output_path": str(resolved_output_path),
            "ctgep_output_path": str(resolved_ctgep_output_path),
            "celltype_mapping": resolved_mapping,
            "samplenum_per_ct": samplenum_per_ct,
            "val_samplenum_per_patient": val_samplenum_per_patient,
            "samplenum_all_train": n_training_samples,
            "samplenum_all_val": n_validation_samples,
            "n_cells": cells_per_sample,
            "uniform_alpha": uniform_alpha,
            "dominant_alphas": [50.0, 15.0, 5.0] if dominant_alphas is None else dominant_alphas,
            "deg_fdr_rate": deg_fdr_rate,
            "deg_num_per_ct": deg_num_per_ct,
            "epochs": epochs,
            "batch_size": batch_size,
            "learning_rate": learning_rate,
            "weight_decay": weight_decay,
            "coeff_ct_gep": coeff_ct_gep,
            "coeff_prop": coeff_prop,
            "sched_step": sched_step,
            "sched_gamma": sched_gamma,
            "preprocess_mode": preprocess_mode,
            "split_mode": split_mode,
            "split_seed": split_seed,
            "seed": seed,
        }
        payload_path.write_text(json.dumps(payload), encoding="utf-8")
        run_command(command_base + [str(runner_path), str(payload_path)], log_path=log_path)
        result = read_sample_by_celltype(str(resolved_output_path))

    output_paths = {}
    if output_path is not None or extra_output_dir is not None:
        output_paths["proportion"] = str(resolved_output_path)
    if ctgep_output_path is not None or extra_output_dir is not None:
        output_paths["ctgep"] = str(resolved_ctgep_output_path)
    if extra_output_dir is not None:
        output_paths["extra_output_dir"] = extra_output_dir

    return DeconvolutionResult(
        tool="blue",
        proportion=result,
        metadata={
            "mixture_transform": mixture_transform,
            "scrna_transform": scrna_transform,
            "python_env": "" if python_env is None else python_env,
            "blue_repo_path": str(resolved_blue_repo),
            "shared_gene_count": len(mixture_aligned.index),
            "input_gene_count": len(resolved_input_genes),
            "output_gene_count": len(resolved_output_genes),
            "cell_type_count": len(resolved_mapping),
            "batch_id_count": int(aligned_batch_ids.nunique()),
            "n_training_samples": n_training_samples,
            "n_validation_samples": n_validation_samples,
            "cells_per_sample": cells_per_sample,
            "epochs": epochs,
            "batch_size": batch_size,
            "learning_rate": learning_rate,
            "preprocess_mode": preprocess_mode,
            "split_mode": split_mode,
            "seed": seed,
            "output_path": "" if output_path is None else output_path,
            "ctgep_output_path": "" if ctgep_output_path is None else ctgep_output_path,
            "extra_output_dir": "" if extra_output_dir is None else extra_output_dir,
            "log_path": "" if log_path is None else log_path,
            "output_semantics": "BLUE predicted proportions; native ctGEP h5ad is recorded in output_paths when persistent",
        },
        output_paths=output_paths,
    )


def run_bisque(
    mixture: pd.DataFrame,
    scrna_counts: pd.DataFrame,
    cell_types: pd.Series,
    batch_ids: pd.Series,
    mixture_transform: str = "raw",
    scrna_transform: str = "raw",
    use_overlap: bool = False,
    old_cpm: bool = True,
) -> DeconvolutionResult:
    """Run Bisque RNA ReferenceBasedDecomposition."""

    validate_expression_matrix(mixture, "mixture")
    validate_expression_matrix(scrna_counts, "scrna_counts")
    _validate_scrna_metadata(scrna_counts, cell_types, batch_ids)

    mixture_norm = normalize_gene_names(mixture)
    scrna_norm = normalize_gene_names(scrna_counts)
    mixture_aligned, scrna_aligned = align_mixture_and_scrna(mixture_norm, scrna_norm)

    mixture_transformed = apply_transform(mixture_aligned, transform=mixture_transform)
    scrna_transformed = apply_transform(scrna_aligned, transform=scrna_transform)

    with tempfile.TemporaryDirectory(prefix="foli_bisque_") as workdir:
        workdir_path = Path(workdir)
        mixture_path = workdir_path / "mixture.tsv"
        scrna_path = workdir_path / "scrna.tsv"
        cell_types_path = workdir_path / "cell_types.tsv"
        batch_ids_path = workdir_path / "batch_ids.tsv"

        write_expression_matrix(mixture_transformed, str(mixture_path))
        write_expression_matrix(scrna_transformed, str(scrna_path))
        write_vector(pd.Series(cell_types.to_numpy()), str(cell_types_path), "cell_type")
        write_vector(pd.Series(batch_ids.to_numpy()), str(batch_ids_path), "batch_id")

        result, _ = run_r_tool(
            tool="bisque",
            payload={
                "mixture_path": str(mixture_path),
                "scrna_path": str(scrna_path),
                "cell_types_path": str(cell_types_path),
                "batch_ids_path": str(batch_ids_path),
                "use_overlap": use_overlap,
                "old_cpm": old_cpm,
            },
        )

    return DeconvolutionResult(
        tool="bisque",
        proportion=result,
        metadata={
            "mixture_transform": mixture_transform,
            "scrna_transform": scrna_transform,
            "use_overlap": use_overlap,
            "old_cpm": old_cpm,
        },
    )


def run_bayesprism(
    mixture: pd.DataFrame,
    scrna_counts: pd.DataFrame,
    cell_types: pd.Series,
    mixture_transform: str = "raw",
    scrna_transform: str = "raw",
    cell_states: pd.Series | None = None,
    tum_key: str | None = None,
    update_gibbs: bool = True,
    n_cores: int = 1,
    which_theta: str = "final",
    state_or_type: str = "type",
    outlier_cut: float = 1.0,
    outlier_fraction: float = 1.0,
    pseudo_min: float = 1e-8,
) -> DeconvolutionResult:
    """Run BayesPrism deconvolution with single-cell references."""

    validate_expression_matrix(mixture, "mixture")
    validate_expression_matrix(scrna_counts, "scrna_counts")
    if len(cell_types) != scrna_counts.shape[1]:
        raise ValueError("cell_types length must match number of scRNA columns")
    if cell_states is not None and len(cell_states) != scrna_counts.shape[1]:
        raise ValueError("cell_states length must match number of scRNA columns")

    mixture_norm = normalize_gene_names(mixture)
    scrna_norm = normalize_gene_names(scrna_counts)
    mixture_aligned, scrna_aligned = align_mixture_and_scrna(mixture_norm, scrna_norm)

    mixture_transformed = apply_transform(mixture_aligned, transform=mixture_transform)
    scrna_transformed = apply_transform(scrna_aligned, transform=scrna_transform)

    with tempfile.TemporaryDirectory(prefix="foli_bayesprism_") as workdir:
        workdir_path = Path(workdir)
        mixture_path = workdir_path / "mixture.tsv"
        scrna_path = workdir_path / "scrna.tsv"
        cell_types_path = workdir_path / "cell_types.tsv"

        write_expression_matrix(mixture_transformed, str(mixture_path))
        write_expression_matrix(scrna_transformed, str(scrna_path))
        write_vector(pd.Series(cell_types.to_numpy()), str(cell_types_path), "cell_type")

        payload: dict[str, object] = {
            "mixture_path": str(mixture_path),
            "scrna_path": str(scrna_path),
            "cell_types_path": str(cell_types_path),
            "tum_key": tum_key,
            "update_gibbs": update_gibbs,
            "n_cores": n_cores,
            "which_theta": which_theta,
            "state_or_type": state_or_type,
            "outlier_cut": outlier_cut,
            "outlier_fraction": outlier_fraction,
            "pseudo_min": pseudo_min,
        }

        if cell_states is not None:
            cell_states_path = workdir_path / "cell_states.tsv"
            write_vector(pd.Series(cell_states.to_numpy()), str(cell_states_path), "cell_state")
            payload["cell_states_path"] = str(cell_states_path)

        result, _ = run_r_tool(tool="bayesprism", payload=payload)

    return DeconvolutionResult(
        tool="bayesprism",
        proportion=result,
        metadata={
            "mixture_transform": mixture_transform,
            "scrna_transform": scrna_transform,
            "update_gibbs": update_gibbs,
            "which_theta": which_theta,
            "state_or_type": state_or_type,
            "outlier_cut": outlier_cut,
            "outlier_fraction": outlier_fraction,
            "pseudo_min": pseudo_min,
        },
    )
