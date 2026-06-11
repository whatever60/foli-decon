"""Run AutoGeneS deconvolution from a JSON payload."""

import json
import sys
from pathlib import Path

import autogenes as ag
import numpy as np
import pandas as pd


def main() -> None:
    """Read a payload, run AutoGeneS feature selection and deconvolution, and write outputs."""

    payload_path = sys.argv[1]
    payload = json.loads(Path(payload_path).read_text(encoding="utf-8"))

    mixture = pd.read_csv(payload["mixture_path"], sep="\t", index_col=0)
    signature = pd.read_csv(payload["signature_path"], sep="\t", index_col=0)

    reference = signature.T
    bulk = mixture.T

    ag.init(reference)
    ag.optimize(
        ngen=int(payload["ngen"]),
        mode=payload["mode"],
        nfeatures=int(payload["nfeatures"]),
        seed=int(payload["seed"]),
        verbose=bool(payload["verbose"]),
    )
    selected_genes = ag.select(index=int(payload["selection_index"]))
    coefficients = ag.deconvolve(bulk, model=payload["model"])

    result = pd.DataFrame(
        coefficients,
        index=bulk.index,
        columns=reference.index,
    )
    if bool(payload["clip_negative"]):
        result[result < 0] = 0.0
    if bool(payload["normalize_coefficients"]):
        row_sums = result.sum(axis=1).to_numpy()
        result = result.div(np.where(row_sums == 0.0, 1.0, row_sums), axis=0)

    output_path = Path(payload["output_path"])
    output_path.parent.mkdir(parents=True, exist_ok=True)
    result.to_csv(output_path, sep="\t")

    selected_features_path = Path(payload["selected_features_path"])
    selected_features_path.parent.mkdir(parents=True, exist_ok=True)
    pd.DataFrame({"gene": pd.Index(selected_genes).astype(str)}).to_csv(
        selected_features_path,
        sep="\t",
        index=False,
    )


if __name__ == "__main__":
    main()
