"""Run TAPE deconvolution from a JSON payload."""

import json
import functools
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")

import pandas as pd
import torch
from TAPE import Deconvolution


def main() -> None:
    """Read a payload, run TAPE, and write sample-by-cell type fractions."""

    payload_path = sys.argv[1]
    payload = json.loads(Path(payload_path).read_text(encoding="utf-8"))

    sc_ref = pd.read_csv(payload["scrna_path"], sep="\t", index_col=0)
    mixture = pd.read_csv(payload["mixture_path"], sep="\t", index_col=0)
    # Native TAPE saves and reloads a model created in this same temp run.
    # PyTorch 2.6 changed torch.load's default to weights_only=True.
    torch.load = functools.partial(torch.load, weights_only=False)
    signature, prediction = Deconvolution(
        sc_ref,
        mixture,
        sep="\t",
        variance_threshold=float(payload["variance_threshold"]),
        scaler=payload["scaler"],
        datatype=payload["datatype"],
        genelenfile=None if payload["gene_lengths_table_path"] == "" else payload["gene_lengths_table_path"],
        d_prior=None,
        mode=payload["mode"],
        adaptive=bool(payload["adaptive"]),
        save_model_name=None,
        sparse=bool(payload["sparse"]),
        batch_size=int(payload["batch_size"]),
        epochs=int(payload["epochs"]),
        seed=int(payload["seed"]),
    )

    output_path = Path(payload["output_path"])
    output_path.parent.mkdir(parents=True, exist_ok=True)
    prediction.to_csv(output_path, sep="\t")

    if payload["signature_output_path"] != "" and hasattr(signature, "to_csv"):
        signature_output_path = Path(payload["signature_output_path"])
        signature_output_path.parent.mkdir(parents=True, exist_ok=True)
        signature.to_csv(signature_output_path, sep="\t")


if __name__ == "__main__":
    main()
