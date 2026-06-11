"""Run the BLUE multi-step deconvolution pipeline from a JSON payload."""

import json
import os
import shutil
import subprocess
import sys
from pathlib import Path

import anndata as ad
import pandas as pd
import scipy.sparse as sp
import yaml


def main() -> None:
    """Stage BLUE inputs, run native pipeline scripts, and collect outputs."""

    payload = json.loads(Path(sys.argv[1]).read_text(encoding="utf-8"))
    blue_repo = Path(payload["blue_repo_path"])
    workdir = Path(payload["workdir"])
    workdir.mkdir(parents=True, exist_ok=True)

    scrna = pd.read_csv(payload["scrna_matrix_path"], sep="\t", index_col=0)
    metadata = pd.read_csv(payload["scrna_metadata_path"], sep="\t")
    metadata.index = metadata["cell_id"].astype(str)

    sc_h5ad_dir = workdir / "single_cell_h5ad"
    sc_h5ad_dir.mkdir(parents=True, exist_ok=True)
    library_ids = sorted(metadata["library_id"].astype(str).unique().tolist())
    for library_index, library_id in enumerate(library_ids):
        library_metadata = metadata.loc[metadata["library_id"].astype(str) == library_id].copy()
        cell_ids = library_metadata.index.astype(str).tolist()
        obs = library_metadata.loc[cell_ids, ["cell_type", "library_id"]].copy()
        obs.index = cell_ids
        adata = ad.AnnData(
            X=sp.csr_matrix(scrna.loc[:, cell_ids].T.to_numpy(dtype=float)),
            obs=obs,
            var=pd.DataFrame(index=scrna.index.astype(str)),
        )
        adata.write_h5ad(sc_h5ad_dir / f"library_{library_index}.h5ad")

    mapping_path = workdir / "celltype_mapping.yaml"
    mapping_path.write_text(yaml.safe_dump(payload["celltype_mapping"], sort_keys=False), encoding="utf-8")

    config = {
        "paths": {
            "sc_h5ad_dir": str(sc_h5ad_dir),
            "sc_h5ad_glob": "*.h5ad",
            "bulk_tsv": payload["mixture_path"],
            "bulk_metadata_tsv": None,
            "combined_sc_h5ad": str(workdir / "combined_labeled.h5ad"),
            "celltype_mapping_yaml": str(mapping_path),
            "gene_lists_dir": str(workdir / "deconv_training"),
            "pseudobulk_dir": str(workdir / "deconv_training" / "pseudobulks"),
            "ckpt_root": str(workdir / "deconv_ckpt"),
            "predictions_dir": str(workdir / "deconv_predictions"),
            "input_gene_list": payload["input_gene_list_path"],
            "output_gene_list": payload["output_gene_list_path"],
        },
        "sc_schema": {
            "fine_celltype_col": "cell_type",
            "library_id_col": "library_id",
            "raw_counts_layer": None,
            "ensembl_var_names": False,
        },
        "pseudobulk": {
            "samplenum_per_ct": payload["samplenum_per_ct"],
            "val_samplenum_per_patient": payload["val_samplenum_per_patient"],
            "samplenum_all_train": payload["samplenum_all_train"],
            "samplenum_all_val": payload["samplenum_all_val"],
            "n_cells": payload["n_cells"],
            "train_ratio": 0.8,
            "val_ratio": 0.1,
            "seed": payload["seed"],
            "uniform_alpha": payload["uniform_alpha"],
            "dominant_alphas": payload["dominant_alphas"],
        },
        "degs": {
            "fdr_rate": payload["deg_fdr_rate"],
            "num_deg_per_ct": payload["deg_num_per_ct"],
        },
        "training": {
            "batch_size": payload["batch_size"],
            "epochs": payload["epochs"],
            "lr": payload["learning_rate"],
            "weight_decay": payload["weight_decay"],
            "coeff_ct_gep": payload["coeff_ct_gep"],
            "coeff_prop": payload["coeff_prop"],
            "sched_step": payload["sched_step"],
            "sched_gamma": payload["sched_gamma"],
            "use_amp": False,
            "seed": payload["seed"],
            "cuda": 0,
            "preprocess_mode": payload["preprocess_mode"],
            "split_mode": payload["split_mode"],
            "split_seed": payload["split_seed"],
        },
        "predict": {
            "cuda": 0,
        },
    }
    config_path = workdir / "pipeline_config.yaml"
    config_path.write_text(yaml.safe_dump(config, sort_keys=False), encoding="utf-8")

    env = os.environ.copy()
    env["CUDA_VISIBLE_DEVICES"] = ""
    command_prefix = [sys.executable]
    script_dir = blue_repo / "scripts" / "deconv"
    commands = [
        command_prefix + [str(script_dir / "00_combine_sc_h5ads.py"), "--config", str(config_path)],
        command_prefix + [str(script_dir / "01_build_gene_lists.py"), "--config", str(config_path)],
        command_prefix + [str(script_dir / "02_simulate_pseudobulk.py"), "--config", str(config_path)],
        command_prefix + [str(script_dir / "03_preprocess_bulk.py"), "--config", str(config_path)],
        command_prefix
        + [
            str(script_dir / "04_train_unet.py"),
            "--config",
            str(config_path),
            "--run-stamp",
            "foli_blue",
        ],
        command_prefix + [str(script_dir / "05_predict_bulk.py"), "--config", str(config_path)],
    ]
    for command in commands:
        process = subprocess.run(command, cwd=blue_repo, env=env, capture_output=True, text=True, check=False)
        print(f"$ {' '.join(command)}")
        print(process.stdout)
        print(process.stderr, file=sys.stderr)
        if process.returncode != 0:
            raise RuntimeError(f"BLUE command failed with exit code {process.returncode}: {' '.join(command)}")

    prediction_dir = workdir / "deconv_predictions" / "foli_blue"
    proportion_candidates = sorted(prediction_dir.glob("predicted_proportions_ep*.csv"))
    ctgep_candidates = sorted(prediction_dir.glob("predicted_ctGEP_ep*.h5ad"))
    native_proportion_path = proportion_candidates[-1]
    native_ctgep_path = ctgep_candidates[-1]

    proportion = pd.read_csv(native_proportion_path, index_col=0)
    Path(payload["output_path"]).parent.mkdir(parents=True, exist_ok=True)
    proportion.to_csv(payload["output_path"], sep="\t")
    Path(payload["ctgep_output_path"]).parent.mkdir(parents=True, exist_ok=True)
    shutil.copyfile(native_ctgep_path, payload["ctgep_output_path"])


if __name__ == "__main__":
    main()
