"""Run MarkerMap feature selection from a JSON payload."""

import json
import sys
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import pytorch_lightning as pl
import torch
from markermap.vae_models import MarkerMap, train_model


def _split_indices(labels: pd.Series, train_fraction: float, seed: int) -> tuple[np.ndarray, np.ndarray]:
    """Create a deterministic stratified train/validation split."""

    rng = np.random.default_rng(seed)
    train_indices = []
    val_indices = []
    for label in sorted(labels.unique()):
        label_indices = np.where(labels.to_numpy() == label)[0]
        rng.shuffle(label_indices)
        n_val = int(np.floor((1.0 - train_fraction) * len(label_indices)))
        if n_val == 0 and len(label_indices) > 2 and train_fraction < 1.0:
            n_val = 1
        val_indices.extend(label_indices[:n_val])
        train_indices.extend(label_indices[n_val:])
    return np.asarray(train_indices, dtype=int), np.asarray(val_indices, dtype=int)


def main() -> None:
    """Read a payload, train MarkerMap, and write selected-gene tables."""

    payload_path = sys.argv[1]
    payload = json.loads(Path(payload_path).read_text(encoding="utf-8"))

    scrna = pd.read_csv(payload["scrna_path"], sep="\t", index_col=0)
    labels_table = pd.read_csv(payload["labels_path"], sep="\t")
    labels = labels_table["cell_type"].astype(str)
    categories = sorted(labels.unique())

    seed = int(payload["seed"])
    np.random.seed(seed)
    torch.manual_seed(seed)
    pl.seed_everything(seed, workers=True, verbose=False)

    adata = ad.AnnData(X=scrna.T.to_numpy(dtype=np.float32))
    adata.obs["cell_type"] = pd.Categorical(labels.to_numpy(), categories=categories)
    adata.var_names = scrna.index.astype(str)
    train_indices, val_indices = _split_indices(labels, float(payload["train_fraction"]), seed)

    train_dataloader, val_dataloader = MarkerMap.prepareData(
        adata,
        train_indices,
        val_indices,
        "cell_type",
        None,
        batch_size=int(payload["batch_size"]),
    )
    model = MarkerMap(
        input_size=adata.shape[1],
        hidden_layer_size=int(payload["hidden_layer_size"]),
        z_size=int(payload["z_size"]),
        num_classes=len(categories),
        k=int(payload["nfeatures"]),
        loss_tradeoff=float(payload["loss_tradeoff"]),
    )
    train_model(
        model,
        train_dataloader,
        val_dataloader,
        min_epochs=int(payload["min_epochs"]),
        max_epochs=int(payload["max_epochs"]),
        auto_lr=bool(payload["auto_lr"]),
        early_stopping_patience=int(payload["early_stopping_patience"]),
        verbose=bool(payload["verbose"]),
    )

    marker_indices = model.markers().clone().cpu().detach().numpy().astype(int)
    selected_genes = scrna.index[marker_indices].astype(str)
    ranking = pd.DataFrame(
        {
            "rank": np.arange(1, len(selected_genes) + 1, dtype=int),
            "gene": selected_genes,
            "feature_index": marker_indices,
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
