"""Run AutoGeneS feature selection from a JSON payload."""

import json
import sys
from pathlib import Path

import autogenes as ag
import pandas as pd


def main() -> None:
    """Read a payload, run AutoGeneS feature selection, and write selected genes."""

    payload_path = sys.argv[1]
    payload = json.loads(Path(payload_path).read_text(encoding="utf-8"))

    signature = pd.read_csv(payload["signature_path"], sep="\t", index_col=0)
    reference = signature.T

    ag.init(reference)
    ag.optimize(
        ngen=int(payload["ngen"]),
        mode=payload["mode"],
        nfeatures=int(payload["nfeatures"]),
        seed=int(payload["seed"]),
        verbose=bool(payload["verbose"]),
    )
    selected_genes = ag.select(index=int(payload["selection_index"]))

    selected_features_path = Path(payload["selected_features_path"])
    selected_features_path.parent.mkdir(parents=True, exist_ok=True)
    pd.DataFrame({"gene": pd.Index(selected_genes).astype(str)}).to_csv(
        selected_features_path,
        sep="\t",
        index=False,
    )


if __name__ == "__main__":
    main()
