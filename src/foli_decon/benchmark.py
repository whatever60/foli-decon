"""Two-split simulation benchmark workflow for deconvolution tools."""

from dataclasses import dataclass
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
from scipy.spatial.distance import jensenshannon
from scipy.sparse import issparse

from foli_decon.api import run_deconvolution


@dataclass
class BenchmarkResult:
    """Container for benchmark outputs."""

    metrics: pd.DataFrame
    truth_by_scenario: dict[str, pd.DataFrame]
    predictions: dict[str, pd.DataFrame]


def _build_truth_scenarios(
    cell_types: list[str],
    n_samples_per_scenario: int,
    dirichlet_alpha: float,
    dominant_fraction: float,
    one_cell_type_fraction: float,
    rng: np.random.Generator,
) -> dict[str, pd.DataFrame]:
    """Create per-sample ground-truth fractions for both simulation scenarios."""

    dirichlet_values = rng.dirichlet(
        np.ones(len(cell_types)) * dirichlet_alpha,
        size=n_samples_per_scenario,
    )
    dirichlet_truth = pd.DataFrame(
        dirichlet_values,
        index=[f"dirichlet_uniform_{i}" for i in range(n_samples_per_scenario)],
        columns=cell_types,
    )

    dominant_values = np.zeros((n_samples_per_scenario, len(cell_types)), dtype=float)
    for i in range(n_samples_per_scenario):
        dominant_index = int(rng.integers(0, len(cell_types)))
        if len(cell_types) == 1:
            dominant_values[i, dominant_index] = 1.0
            continue
        non_dominant = [j for j in range(len(cell_types)) if j != dominant_index]
        remaining = rng.dirichlet(np.ones(len(non_dominant))) * (1.0 - dominant_fraction)
        dominant_values[i, dominant_index] = dominant_fraction
        dominant_values[i, non_dominant] = remaining

    dominant_truth = pd.DataFrame(
        dominant_values,
        index=[f"one_type_dominant_{i}" for i in range(n_samples_per_scenario)],
        columns=cell_types,
    )

    one_cell_type_values = np.zeros((n_samples_per_scenario, len(cell_types)), dtype=float)
    for i in range(n_samples_per_scenario):
        dominant_index = i % len(cell_types)
        if len(cell_types) == 1:
            one_cell_type_values[i, dominant_index] = 1.0
            continue
        if one_cell_type_fraction >= 1.0:
            one_cell_type_values[i, dominant_index] = 1.0
            continue
        non_dominant = [j for j in range(len(cell_types)) if j != dominant_index]
        remaining = rng.dirichlet(np.ones(len(non_dominant))) * (1.0 - one_cell_type_fraction)
        one_cell_type_values[i, dominant_index] = one_cell_type_fraction
        one_cell_type_values[i, non_dominant] = remaining

    one_cell_type_truth = pd.DataFrame(
        one_cell_type_values,
        index=[f"one_cell_type_dominate_{i}" for i in range(n_samples_per_scenario)],
        columns=cell_types,
    )

    return {
        "dirichlet_uniform": dirichlet_truth,
        "one_type_dominant": dominant_truth,
        "one_cell_type_dominate": one_cell_type_truth,
    }


def _simulate_pseudobulk(
    evaluation_counts: np.ndarray,
    evaluation_cell_types: np.ndarray,
    gene_names: pd.Index,
    truth: pd.DataFrame,
    cells_per_bulk: int,
    rng: np.random.Generator,
) -> pd.DataFrame:
    """Simulate gene-by-sample pseudobulk counts from evaluation split cells."""

    cell_types = truth.columns.tolist()
    indices_by_type: dict[str, np.ndarray] = {}
    for cell_type in cell_types:
        indices_by_type[cell_type] = np.flatnonzero(evaluation_cell_types == cell_type)

    bulk_rows = []
    for sample_name in truth.index:
        fractions = truth.loc[sample_name, cell_types].to_numpy(dtype=float)
        sampled_counts = rng.multinomial(cells_per_bulk, fractions)
        sampled_indices = []
        for i, cell_type in enumerate(cell_types):
            n_cells = int(sampled_counts[i])
            if n_cells == 0:
                continue
            pool = indices_by_type[cell_type]
            sampled = rng.choice(pool, size=n_cells, replace=True)
            sampled_indices.append(sampled)
        combined_indices = np.concatenate(sampled_indices)
        if issparse(evaluation_counts):
            bulk_profile = np.asarray(evaluation_counts[combined_indices].sum(axis=0)).ravel()
        else:
            bulk_profile = evaluation_counts[combined_indices].sum(axis=0)
        bulk_rows.append(bulk_profile)

    bulk_matrix = np.vstack(bulk_rows)
    return pd.DataFrame(bulk_matrix.T, index=gene_names, columns=truth.index)


def _build_signature_matrix(
    reference_counts: np.ndarray,
    reference_cell_types: np.ndarray,
    cell_types: list[str],
    gene_names: pd.Index,
) -> pd.DataFrame:
    """Average reference split cells into a gene-by-cell-type signature matrix."""

    columns = []
    for cell_type in cell_types:
        indices = np.flatnonzero(reference_cell_types == cell_type)
        if issparse(reference_counts):
            profile = np.asarray(reference_counts[indices].mean(axis=0)).ravel()
        else:
            profile = reference_counts[indices].mean(axis=0)
        columns.append(profile)
    signature = np.vstack(columns).T
    return pd.DataFrame(signature, index=gene_names, columns=cell_types)


def _to_gene_by_cell_dataframe(
    counts: np.ndarray,
    gene_names: pd.Index,
    cell_names: list[str],
) -> pd.DataFrame:
    """Convert cell-by-gene count matrix into gene-by-cell pandas DataFrame."""

    if issparse(counts):
        values = counts.toarray().T
    else:
        values = np.asarray(counts).T
    return pd.DataFrame(values, index=gene_names, columns=cell_names)


def _subset_panel_genes(
    gene_by_sample: pd.DataFrame,
    panel_genes: list[str] | None,
) -> pd.DataFrame:
    """Subset a gene-indexed matrix to panel genes while preserving original order."""

    if panel_genes is None:
        return gene_by_sample
    panel = set(pd.Index(panel_genes).astype(str))
    shared = [gene for gene in gene_by_sample.index if gene in panel]
    return gene_by_sample.loc[shared]


def _collapse_duplicated_genes(gene_by_sample: pd.DataFrame) -> pd.DataFrame:
    """Collapse duplicated gene rows by summation."""

    if gene_by_sample.index.is_unique:
        return gene_by_sample
    return gene_by_sample.groupby(level=0, sort=False).sum()


def _normalize_and_collapse_gene_names(gene_by_sample: pd.DataFrame) -> pd.DataFrame:
    """Replace dots with dashes in gene names and collapse duplicates."""

    normalized = gene_by_sample.copy()
    normalized.index = normalized.index.astype(str).str.replace(".", "-", regex=False)
    return _collapse_duplicated_genes(normalized)


def _normalize_rows(matrix: np.ndarray) -> np.ndarray:
    """Normalize rows to sum to one."""

    row_sums = matrix.sum(axis=1, keepdims=True)
    if (row_sums <= 0).any():
        raise ValueError("Each sample row must have positive total mass for normalization")
    return matrix / row_sums


def _compute_metrics(
    truth: pd.DataFrame,
    prediction: pd.DataFrame,
) -> dict[str, float | int]:
    """Compute accuracy metrics after sample/cell-type alignment."""

    shared_samples = [sample for sample in truth.index if sample in prediction.index]
    shared_cell_types = [cell_type for cell_type in truth.columns if cell_type in prediction.columns]
    if len(shared_samples) == 0:
        raise ValueError("No overlapping sample names between truth and prediction")
    if len(shared_cell_types) == 0:
        raise ValueError("No overlapping cell-type columns between truth and prediction")

    truth_values = truth.loc[shared_samples, shared_cell_types].to_numpy(dtype=float)
    pred_values = prediction.loc[shared_samples, shared_cell_types].to_numpy(dtype=float)
    if not np.isfinite(truth_values).all():
        raise ValueError("Truth matrix contains non-finite values")
    if not np.isfinite(pred_values).all():
        raise ValueError("Prediction matrix contains non-finite values")
    truth_values = np.clip(truth_values, 0.0, None)
    pred_values = np.clip(pred_values, 0.0, None)
    truth_norm = _normalize_rows(truth_values)
    pred_norm = _normalize_rows(pred_values)
    if not np.isfinite(truth_norm).all():
        raise ValueError("Normalized truth matrix contains non-finite values")
    if not np.isfinite(pred_norm).all():
        raise ValueError("Normalized prediction matrix contains non-finite values")

    diff = pred_norm - truth_norm
    mae = float(np.abs(diff).mean())
    rmse = float(np.sqrt((diff**2).mean()))

    js_values = []
    cosine_values = []
    kl_truth_to_pred_values = []
    kl_pred_to_truth_values = []
    sym_kl_values = []
    hellinger_values = []
    total_variation_values = []
    bray_curtis_values = []
    epsilon = 1e-12
    for i in range(truth_norm.shape[0]):
        truth_i = np.clip(truth_norm[i], epsilon, None)
        pred_i = np.clip(pred_norm[i], epsilon, None)
        truth_i = truth_i / truth_i.sum()
        pred_i = pred_i / pred_i.sum()

        js_values.append(float(jensenshannon(truth_i, pred_i) ** 2))
        norm_truth = float(np.linalg.norm(truth_i))
        norm_pred = float(np.linalg.norm(pred_i))
        cosine_values.append(float(np.dot(truth_i, pred_i) / (norm_truth * norm_pred)))
        kl_truth_to_pred = float(np.sum(truth_i * np.log(truth_i / pred_i)))
        kl_pred_to_truth = float(np.sum(pred_i * np.log(pred_i / truth_i)))
        kl_truth_to_pred_values.append(kl_truth_to_pred)
        kl_pred_to_truth_values.append(kl_pred_to_truth)
        sym_kl_values.append(float(0.5 * (kl_truth_to_pred + kl_pred_to_truth)))
        hellinger_values.append(
            float(np.sqrt(0.5 * np.sum((np.sqrt(truth_i) - np.sqrt(pred_i)) ** 2)))
        )
        total_variation_values.append(float(0.5 * np.sum(np.abs(truth_i - pred_i))))
        bray_curtis_values.append(float(np.sum(np.abs(truth_i - pred_i)) / np.sum(truth_i + pred_i)))

    return {
        "mae": mae,
        "rmse": rmse,
        "js_divergence": float(np.mean(js_values)),
        "cosine_similarity": float(np.mean(cosine_values)),
        "kl_divergence_truth_to_pred": float(np.mean(kl_truth_to_pred_values)),
        "kl_divergence_pred_to_truth": float(np.mean(kl_pred_to_truth_values)),
        "kl_divergence_symmetric": float(np.mean(sym_kl_values)),
        "hellinger_distance": float(np.mean(hellinger_values)),
        "total_variation_distance": float(np.mean(total_variation_values)),
        "bray_curtis_distance": float(np.mean(bray_curtis_values)),
        "n_shared_samples": len(shared_samples),
        "n_shared_cell_types": len(shared_cell_types),
    }


def _prepare_split_indices(
    cell_types: np.ndarray,
    selected_cell_types: list[str],
    reference_fraction: float,
    max_reference_cells_per_type: int,
    max_evaluation_cells_per_type: int,
    rng: np.random.Generator,
    split_strategy: str,
    split_keys: np.ndarray | None,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Create reference/evaluation split indices stratified by cell type."""

    if split_strategy not in {"random_cell_type", "key_stratified"}:
        raise ValueError("split_strategy must be one of {'random_cell_type', 'key_stratified'}")
    if split_strategy == "key_stratified" and split_keys is None:
        raise ValueError("split_keys is required when split_strategy='key_stratified'")

    reference_key_set: set[str] = set()
    if split_strategy == "key_stratified":
        selected_mask = np.isin(cell_types, selected_cell_types)
        selected_indices = np.flatnonzero(selected_mask)
        unique_keys = np.unique(split_keys[selected_indices])
        if len(unique_keys) < 2:
            raise ValueError("key_stratified split requires at least two unique key values")
        shuffled_keys = rng.permutation(unique_keys)
        n_reference_keys = int(round(len(unique_keys) * reference_fraction))
        if n_reference_keys < 1:
            n_reference_keys = 1
        if n_reference_keys >= len(unique_keys):
            n_reference_keys = len(unique_keys) - 1
        reference_key_set = set(shuffled_keys[:n_reference_keys].tolist())

    reference_indices_parts = []
    evaluation_indices_parts = []
    reference_labels_parts = []
    evaluation_labels_parts = []

    for cell_type in selected_cell_types:
        indices = np.flatnonzero(cell_types == cell_type)
        if split_strategy == "random_cell_type":
            shuffled = rng.permutation(indices)
            n_reference = int(round(len(indices) * reference_fraction))
            if n_reference < 1:
                n_reference = 1
            if n_reference >= len(indices):
                n_reference = len(indices) - 1
            if n_reference < 1:
                raise ValueError(f"Cell type '{cell_type}' does not have enough cells for two-way split")
            reference_indices = shuffled[:n_reference]
            evaluation_indices = shuffled[n_reference:]
        else:
            reference_mask = np.array([split_keys[index] in reference_key_set for index in indices])
            reference_indices = indices[reference_mask]
            evaluation_indices = indices[~reference_mask]
            if len(reference_indices) < 1 or len(evaluation_indices) < 1:
                raise ValueError(
                    f"Cell type '{cell_type}' has empty split under key_stratified strategy"
                )

        if len(reference_indices) > max_reference_cells_per_type:
            reference_indices = rng.choice(reference_indices, size=max_reference_cells_per_type, replace=False)
        if len(evaluation_indices) > max_evaluation_cells_per_type:
            evaluation_indices = rng.choice(
                evaluation_indices,
                size=max_evaluation_cells_per_type,
                replace=False,
            )

        reference_indices_parts.append(reference_indices)
        evaluation_indices_parts.append(evaluation_indices)
        reference_labels_parts.append(np.repeat(cell_type, len(reference_indices)))
        evaluation_labels_parts.append(np.repeat(cell_type, len(evaluation_indices)))

    reference_indices = np.concatenate(reference_indices_parts)
    evaluation_indices = np.concatenate(evaluation_indices_parts)
    reference_labels = np.concatenate(reference_labels_parts)
    evaluation_labels = np.concatenate(evaluation_labels_parts)
    return reference_indices, evaluation_indices, reference_labels, evaluation_labels


def run_two_split_benchmark_from_adata(
    adata: ad.AnnData,
    dataset_name: str,
    tools: list[str],
    cell_type_col: str,
    batch_col: str,
    count_layer: str | None = "raw_counts",
    reference_fraction: float = 0.5,
    min_cells_per_type: int = 30,
    max_cell_types: int | None = None,
    max_reference_cells_per_type: int = 200,
    max_evaluation_cells_per_type: int = 500,
    n_samples_per_scenario: int = 24,
    cells_per_bulk: int = 400,
    dirichlet_alpha: float = 1.0,
    dominant_fraction: float = 0.8,
    one_cell_type_fraction: float = 1.0,
    split_strategy: str = "key_stratified",
    split_key_col: str | None = None,
    panel_genes: list[str] | None = None,
    xcell2_object_path: str | None = None,
    cibersortx_username: str | None = None,
    cibersortx_token: str | None = None,
    tool_kwargs: dict[str, dict[str, object]] | None = None,
    seed: int = 0,
) -> BenchmarkResult:
    """Run two-split pseudobulk benchmark against selected tools."""

    rng = np.random.default_rng(seed)
    cell_types = adata.obs[cell_type_col].astype(str).to_numpy()
    batches = adata.obs[batch_col].astype(str).to_numpy()
    split_keys = None
    if split_strategy == "key_stratified":
        key_col = split_key_col
        if key_col is None:
            key_col = batch_col
        split_keys = adata.obs[key_col].astype(str).to_numpy()

    cell_type_counts = pd.Series(cell_types).value_counts()
    cell_type_counts = cell_type_counts[cell_type_counts >= min_cells_per_type]
    if max_cell_types is not None:
        cell_type_counts = cell_type_counts.iloc[:max_cell_types]
    selected_cell_types = cell_type_counts.index.tolist()
    if len(selected_cell_types) == 0:
        raise ValueError("No cell types satisfy min_cells_per_type")

    reference_indices, evaluation_indices, reference_labels, evaluation_labels = _prepare_split_indices(
        cell_types=cell_types,
        selected_cell_types=selected_cell_types,
        reference_fraction=reference_fraction,
        max_reference_cells_per_type=max_reference_cells_per_type,
        max_evaluation_cells_per_type=max_evaluation_cells_per_type,
        rng=rng,
        split_strategy=split_strategy,
        split_keys=split_keys,
    )

    if count_layer is None:
        reference_counts = adata.X[reference_indices]
        evaluation_counts = adata.X[evaluation_indices]
    else:
        reference_counts = adata.layers[count_layer][reference_indices]
        evaluation_counts = adata.layers[count_layer][evaluation_indices]

    gene_names = adata.var_names.astype(str)
    signature = _build_signature_matrix(
        reference_counts=reference_counts,
        reference_cell_types=reference_labels,
        cell_types=selected_cell_types,
        gene_names=gene_names,
    )

    reference_cell_names = [f"ref_cell_{i}" for i in range(len(reference_indices))]
    scrna_counts = _to_gene_by_cell_dataframe(
        counts=reference_counts,
        gene_names=gene_names,
        cell_names=reference_cell_names,
    )
    signature = _normalize_and_collapse_gene_names(signature)
    scrna_counts = _normalize_and_collapse_gene_names(scrna_counts)
    signature = _subset_panel_genes(signature, panel_genes)
    scrna_counts = _subset_panel_genes(scrna_counts, panel_genes)

    scenario_truth = _build_truth_scenarios(
        cell_types=selected_cell_types,
        n_samples_per_scenario=n_samples_per_scenario,
        dirichlet_alpha=dirichlet_alpha,
        dominant_fraction=dominant_fraction,
        one_cell_type_fraction=one_cell_type_fraction,
        rng=rng,
    )

    scenario_mixtures: dict[str, pd.DataFrame] = {}
    for scenario_name, truth in scenario_truth.items():
        mixture = _simulate_pseudobulk(
            evaluation_counts=evaluation_counts,
            evaluation_cell_types=evaluation_labels,
            gene_names=gene_names,
            truth=truth,
            cells_per_bulk=cells_per_bulk,
            rng=rng,
        )
        mixture = _normalize_and_collapse_gene_names(mixture)
        scenario_mixtures[scenario_name] = _subset_panel_genes(mixture, panel_genes)

    reference_cell_types_series = pd.Series(reference_labels, index=reference_cell_names, name="cell_type")
    reference_batches = pd.Series(batches[reference_indices], index=reference_cell_names, name="batch_id")

    rows = []
    predictions: dict[str, pd.DataFrame] = {}
    for tool in tools:
        for scenario_name, mixture in scenario_mixtures.items():
            kwargs: dict[str, object] = {"mixture": mixture}
            if tool in {"epic", "dtangle", "deconrnaseq", "dwls"}:
                kwargs["signature"] = signature
            if tool in {"music", "bisque", "bayesprism"}:
                kwargs["scrna_counts"] = scrna_counts
                kwargs["cell_types"] = reference_cell_types_series
            if tool in {"music", "bisque"}:
                kwargs["batch_ids"] = reference_batches
            if tool == "xcell2":
                if xcell2_object_path is None:
                    raise ValueError("xcell2_object_path is required for tool='xcell2'")
                kwargs["xcell2_object_path"] = xcell2_object_path
            if tool == "cibersortx":
                if cibersortx_username is None or cibersortx_token is None:
                    raise ValueError("cibersortx_username and cibersortx_token are required for tool='cibersortx'")
                kwargs["signature"] = signature
                kwargs["username"] = cibersortx_username
                kwargs["token"] = cibersortx_token

            if tool_kwargs is not None and tool in tool_kwargs:
                kwargs.update(tool_kwargs[tool])

            try:
                result = run_deconvolution(tool=tool, **kwargs)
                predictions[f"{tool}:{scenario_name}"] = result.proportions
                if tool == "xcell2":
                    shared_samples = [
                        sample
                        for sample in scenario_truth[scenario_name].index
                        if sample in result.proportions.index
                    ]
                    shared_cell_types = [
                        cell_type
                        for cell_type in scenario_truth[scenario_name].columns
                        if cell_type in result.proportions.columns
                    ]
                    rows.append(
                        {
                            "dataset": dataset_name,
                            "tool": tool,
                            "scenario": scenario_name,
                            "status": "ok",
                            "error_message": "score_output_not_proportion_metrics_skipped",
                            "mae": np.nan,
                            "rmse": np.nan,
                            "js_divergence": np.nan,
                            "cosine_similarity": np.nan,
                            "kl_divergence_truth_to_pred": np.nan,
                            "kl_divergence_pred_to_truth": np.nan,
                            "kl_divergence_symmetric": np.nan,
                            "hellinger_distance": np.nan,
                            "total_variation_distance": np.nan,
                            "bray_curtis_distance": np.nan,
                            "n_shared_samples": len(shared_samples),
                            "n_shared_cell_types": len(shared_cell_types),
                        }
                    )
                else:
                    metric_values = _compute_metrics(
                        truth=scenario_truth[scenario_name],
                        prediction=result.proportions,
                    )
                    rows.append(
                        {
                            "dataset": dataset_name,
                            "tool": tool,
                            "scenario": scenario_name,
                            "status": "ok",
                            "error_message": "",
                            **metric_values,
                        }
                    )
            except Exception as error:
                rows.append(
                    {
                        "dataset": dataset_name,
                        "tool": tool,
                        "scenario": scenario_name,
                        "status": "failed",
                        "error_message": str(error),
                        "mae": np.nan,
                        "rmse": np.nan,
                        "js_divergence": np.nan,
                        "cosine_similarity": np.nan,
                        "kl_divergence_truth_to_pred": np.nan,
                        "kl_divergence_pred_to_truth": np.nan,
                        "kl_divergence_symmetric": np.nan,
                        "hellinger_distance": np.nan,
                        "total_variation_distance": np.nan,
                        "bray_curtis_distance": np.nan,
                        "n_shared_samples": 0,
                        "n_shared_cell_types": 0,
                    }
                )

    metrics = pd.DataFrame(rows)
    return BenchmarkResult(
        metrics=metrics,
        truth_by_scenario=scenario_truth,
        predictions=predictions,
    )


def run_two_split_benchmark_from_h5ad(
    adata_path: str,
    tools: list[str],
    cell_type_col: str,
    batch_col: str,
    count_layer: str | None = "raw_counts",
    reference_fraction: float = 0.5,
    min_cells_per_type: int = 30,
    max_cell_types: int | None = None,
    max_reference_cells_per_type: int = 200,
    max_evaluation_cells_per_type: int = 500,
    n_samples_per_scenario: int = 24,
    cells_per_bulk: int = 400,
    dirichlet_alpha: float = 1.0,
    dominant_fraction: float = 0.8,
    one_cell_type_fraction: float = 1.0,
    split_strategy: str = "key_stratified",
    split_key_col: str | None = None,
    panel_genes: list[str] | None = None,
    xcell2_object_path: str | None = None,
    cibersortx_username: str | None = None,
    cibersortx_token: str | None = None,
    tool_kwargs: dict[str, dict[str, object]] | None = None,
    seed: int = 0,
) -> BenchmarkResult:
    """Load a dataset from disk, subset cells, and run two-split benchmark."""

    rng = np.random.default_rng(seed)
    dataset_name = Path(adata_path).stem
    backed = ad.read_h5ad(adata_path, backed="r")
    observed_cell_types = backed.obs[cell_type_col].astype(str)
    cell_type_counts = observed_cell_types.value_counts()
    cell_type_counts = cell_type_counts[cell_type_counts >= min_cells_per_type]
    if max_cell_types is not None:
        cell_type_counts = cell_type_counts.iloc[:max_cell_types]
    selected_cell_types = cell_type_counts.index.tolist()
    if len(selected_cell_types) == 0:
        raise ValueError("No cell types satisfy min_cells_per_type in provided h5ad")

    max_cells_per_type = max_reference_cells_per_type + max_evaluation_cells_per_type
    selected_indices = []
    for cell_type in selected_cell_types:
        indices = np.flatnonzero(observed_cell_types.to_numpy() == cell_type)
        if len(indices) > max_cells_per_type:
            indices = rng.choice(indices, size=max_cells_per_type, replace=False)
        selected_indices.append(indices)
    selected_indices_array = np.concatenate(selected_indices)

    adata = backed[selected_indices_array, :].to_memory()
    backed.file.close()
    return run_two_split_benchmark_from_adata(
        adata=adata,
        dataset_name=dataset_name,
        tools=tools,
        cell_type_col=cell_type_col,
        batch_col=batch_col,
        count_layer=count_layer,
        reference_fraction=reference_fraction,
        min_cells_per_type=min_cells_per_type,
        max_cell_types=max_cell_types,
        max_reference_cells_per_type=max_reference_cells_per_type,
        max_evaluation_cells_per_type=max_evaluation_cells_per_type,
        n_samples_per_scenario=n_samples_per_scenario,
        cells_per_bulk=cells_per_bulk,
        dirichlet_alpha=dirichlet_alpha,
        dominant_fraction=dominant_fraction,
        one_cell_type_fraction=one_cell_type_fraction,
        split_strategy=split_strategy,
        split_key_col=split_key_col,
        panel_genes=panel_genes,
        xcell2_object_path=xcell2_object_path,
        cibersortx_username=cibersortx_username,
        cibersortx_token=cibersortx_token,
        tool_kwargs=tool_kwargs,
        seed=seed,
    )
