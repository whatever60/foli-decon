"""Run two-split deconvolution benchmarks for one or more scRNA-seq references."""

import argparse
from pathlib import Path

import pandas as pd

from foli_decon import run_two_split_benchmark_from_h5ad, supported_tools


def _parse_tool_list(raw: str) -> list[str]:
    """Parse comma-separated tool names."""

    return [value.strip() for value in raw.split(",") if value.strip()]


def _load_panel_genes(panel_genes_path: str | None) -> list[str] | None:
    """Load panel genes from first column of a TSV/CSV file."""

    if panel_genes_path is None:
        return None
    panel_df = pd.read_csv(panel_genes_path, sep=r"[\t,]", engine="python")
    return panel_df.iloc[:, 0].astype(str).tolist()


def main() -> None:
    """Parse arguments and execute benchmark pipeline."""

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--adata-paths",
        nargs="+",
        required=True,
        help="One or more .h5ad paths",
    )
    parser.add_argument(
        "--tools",
        default=",".join(supported_tools()),
        help="Comma-separated tool list",
    )
    parser.add_argument("--cell-type-col", default="celltype")
    parser.add_argument("--batch-col", default="donor")
    parser.add_argument("--count-layer", default="raw_counts")
    parser.add_argument("--reference-fraction", type=float, default=0.5)
    parser.add_argument("--min-cells-per-type", type=int, default=30)
    parser.add_argument("--max-cell-types", type=int, default=20)
    parser.add_argument("--max-reference-cells-per-type", type=int, default=200)
    parser.add_argument("--max-evaluation-cells-per-type", type=int, default=500)
    parser.add_argument("--n-samples-per-scenario", type=int, default=24)
    parser.add_argument("--cells-per-bulk", type=int, default=400)
    parser.add_argument("--dirichlet-alpha", type=float, default=1.0)
    parser.add_argument("--dominant-fraction", type=float, default=0.8)
    parser.add_argument("--one-cell-type-fraction", type=float, default=1.0)
    parser.add_argument(
        "--split-strategy",
        choices=["random_cell_type", "key_stratified"],
        default="key_stratified",
    )
    parser.add_argument("--split-key-col")
    parser.add_argument("--panel-genes-path")
    parser.add_argument("--xcell2-object-path")
    parser.add_argument("--xcell2-min-shared-genes", type=float, default=0.2)
    parser.add_argument("--cibersortx-username")
    parser.add_argument("--cibersortx-token")
    parser.add_argument("--seed", type=int, default=0)
    parser.add_argument("--output-dir", required=True)
    args = parser.parse_args()

    tools = _parse_tool_list(args.tools)
    panel_genes = _load_panel_genes(args.panel_genes_path)
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    tool_kwargs: dict[str, dict[str, object]] = {}
    tool_kwargs["xcell2"] = {"min_shared_genes": args.xcell2_min_shared_genes}

    all_metrics = []
    for adata_path in args.adata_paths:
        result = run_two_split_benchmark_from_h5ad(
            adata_path=adata_path,
            tools=tools,
            cell_type_col=args.cell_type_col,
            batch_col=args.batch_col,
            count_layer=args.count_layer,
            reference_fraction=args.reference_fraction,
            min_cells_per_type=args.min_cells_per_type,
            max_cell_types=args.max_cell_types,
            max_reference_cells_per_type=args.max_reference_cells_per_type,
            max_evaluation_cells_per_type=args.max_evaluation_cells_per_type,
            n_samples_per_scenario=args.n_samples_per_scenario,
            cells_per_bulk=args.cells_per_bulk,
            dirichlet_alpha=args.dirichlet_alpha,
            dominant_fraction=args.dominant_fraction,
            one_cell_type_fraction=args.one_cell_type_fraction,
            split_strategy=args.split_strategy,
            split_key_col=args.split_key_col,
            panel_genes=panel_genes,
            xcell2_object_path=args.xcell2_object_path,
            cibersortx_username=args.cibersortx_username,
            cibersortx_token=args.cibersortx_token,
            tool_kwargs=tool_kwargs,
            seed=args.seed,
        )
        dataset_name = Path(adata_path).stem
        metrics_path = output_dir / f"{dataset_name}_metrics.tsv"
        result.metrics.to_csv(metrics_path, sep="\t", index=False)
        all_metrics.append(result.metrics)

    combined = pd.concat(all_metrics, axis=0, ignore_index=True)
    combined.to_csv(output_dir / "combined_metrics.tsv", sep="\t", index=False)


if __name__ == "__main__":
    main()
