"""Run scGeneFit feature selection from a JSON payload."""

import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scGeneFit.functions import get_markers


def main() -> None:
    """Read a payload, run scGeneFit, and write selected-gene tables."""

    payload_path = sys.argv[1]
    payload = json.loads(Path(payload_path).read_text(encoding="utf-8"))

    scrna = pd.read_csv(payload["scrna_path"], sep="\t", index_col=0)
    labels_table = pd.read_csv(payload["labels_path"], sep="\t")
    labels = labels_table["cell_type"].astype(str)
    categories = sorted(labels.unique())
    label_codes = pd.Categorical(labels, categories=categories).codes

    np.random.seed(int(payload["seed"]))
    marker_indices = get_markers(
        scrna.T.to_numpy(dtype=float),
        label_codes,
        int(payload["nfeatures"]),
        method=payload["method"],
        epsilon=float(payload["epsilon"]),
        sampling_rate=float(payload["sampling_rate"]),
        n_neighbors=int(payload["n_neighbors"]),
        max_constraints=int(payload["max_constraints"]),
        redundancy=float(payload["redundancy"]),
        verbose=bool(payload["verbose"]),
    )
    selected_genes = scrna.index[np.asarray(marker_indices, dtype=int)].astype(str)
    ranking = pd.DataFrame(
        {
            "rank": np.arange(1, len(selected_genes) + 1, dtype=int),
            "gene": selected_genes,
            "feature_index": np.asarray(marker_indices, dtype=int),
        }
    )

    selected_features_path = Path(payload["selected_features_path"])
    selected_features_path.parent.mkdir(parents=True, exist_ok=True)
    ranking[["gene"]].to_csv(selected_features_path, sep="\t", index=False)

    ranking_path = Path(payload["ranking_path"])
    ranking_path.parent.mkdir(parents=True, exist_ok=True)
    ranking.to_csv(ranking_path, sep="\t", index=False)


if __name__ == "__main__":
    main()
