"""Run DISSECT fraction deconvolution from a JSON payload."""

import copy
import json
import sys
from pathlib import Path

import anndata as ad
import pandas as pd
import scipy.sparse as sp
from dissect.PropsSimulator.simulator import simulate
from dissect.configs.config import config
from dissect.dissect_frac import run_dissect_frac
from dissect.prepare_data import dataset


def build_config(payload: dict[str, object]) -> dict[str, object]:
    """Translate a foli-decon payload into a DISSECT config dictionary."""

    run_config = copy.deepcopy(config)
    run_config["experiment_folder"] = payload["experiment_folder"]

    run_config["simulation_params"]["scdata"] = payload["scrna_h5ad_path"]
    run_config["simulation_params"]["save_expr"] = False
    run_config["simulation_params"]["n_samples"] = payload["n_training_samples"]
    run_config["simulation_params"]["type"] = "bulk"
    run_config["simulation_params"]["celltype_col"] = "cell_type"
    run_config["simulation_params"]["batch_col"] = payload["batch_col"]
    run_config["simulation_params"]["cells_per_sample"] = payload["cells_per_sample"]
    run_config["simulation_params"]["preprocess"] = None
    run_config["simulation_params"]["filter"]["min_genes"] = payload["min_genes"]
    run_config["simulation_params"]["filter"]["min_cells"] = payload["min_cells"]
    run_config["simulation_params"]["filter"]["mt_cutoff"] = payload["mt_cutoff"]
    run_config["simulation_params"]["filter"]["min_expr"] = payload["min_expr"]
    run_config["simulation_params"]["concentration"] = None
    run_config["simulation_params"]["prop_sparse"] = payload["prop_sparse"]
    run_config["simulation_params"]["generate_component_figures"] = False

    run_config["deconv_params"]["test_dataset"] = payload["mixture_path"]
    run_config["deconv_params"]["test_dataset_format"] = "txt"
    run_config["deconv_params"]["test_dataset_type"] = "bulk"
    run_config["deconv_params"]["duplicated"] = "first"
    run_config["deconv_params"]["normalize_simulated"] = payload["normalize_simulated"]
    run_config["deconv_params"]["normalize_test"] = payload["normalize_test"]
    run_config["deconv_params"]["var_cutoff"] = payload["var_cutoff"]
    run_config["deconv_params"]["test_in_mix"] = payload["test_in_mix"]
    run_config["deconv_params"]["simulated"] = True
    run_config["deconv_params"]["sig_matrix"] = False
    run_config["deconv_params"]["mix"] = payload["mix"]
    run_config["deconv_params"]["save_config"] = True
    run_config["deconv_params"]["network_params"]["n_steps"] = payload["train_steps"]
    run_config["deconv_params"]["network_params"]["lr"] = payload["learning_rate"]
    run_config["deconv_params"]["network_params"]["batch_size"] = payload["batch_size"]
    run_config["deconv_params"]["network_params"]["n_steps_expr"] = payload["train_steps"]
    run_config["deconv_params"]["alpha_range"] = payload["alpha_range"]
    run_config["deconv_params"]["normalization_per_batch"] = payload["normalization_per_batch"]
    run_config["deconv_params"]["models"] = list(range(int(payload["n_models"])))
    return run_config


def write_scrna_h5ad(payload: dict[str, object]) -> None:
    """Write a DISSECT-readable h5ad using the sidecar AnnData version."""

    scrna = pd.read_csv(payload["scrna_matrix_path"], sep="\t", index_col=0)
    metadata = pd.read_csv(payload["scrna_metadata_path"], sep="\t")
    metadata.index = metadata["cell_id"].astype(str)
    adata = ad.AnnData(
        X=sp.csr_matrix(scrna.T.to_numpy(dtype=float)),
        obs=metadata,
        var=pd.DataFrame(index=scrna.index.astype(str)),
    )
    output_path = Path(payload["scrna_h5ad_path"])
    output_path.parent.mkdir(parents=True, exist_ok=True)
    adata.write_h5ad(output_path)


def main() -> None:
    """Read a payload, run DISSECT simulation/preprocessing/fraction estimation, and exit."""

    payload_path = sys.argv[1]
    payload = json.loads(open(payload_path, encoding="utf-8").read())
    write_scrna_h5ad(payload)
    run_config = build_config(payload)
    simulate(run_config)
    dataset(run_config)
    run_dissect_frac(run_config)


if __name__ == "__main__":
    main()
