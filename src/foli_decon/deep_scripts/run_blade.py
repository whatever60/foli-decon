"""Run BLADE deconvolution from a JSON payload."""

import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from Deconvolution.BLADE import Framework


def _as_numeric_list(payload: dict[str, object], key: str) -> list[float]:
    """Return a payload list as floats."""

    return [float(value) for value in payload[key]]


def main() -> None:
    """Read a payload, run BLADE, and write sample-by-cell type fractions."""

    payload_path = sys.argv[1]
    payload = json.loads(Path(payload_path).read_text(encoding="utf-8"))

    mean = pd.read_csv(payload["signature_mean_path"], sep="\t", index_col=0)
    sd = pd.read_csv(payload["signature_sd_path"], sep="\t", index_col=0)
    mixture = pd.read_csv(payload["mixture_path"], sep="\t", index_col=0)

    np.random.seed(int(payload["seed"]))
    final_obj, best_obj, best_set, all_out = Framework(
        mean.to_numpy(),
        sd.to_numpy(),
        mixture.to_numpy(),
        Alphas=_as_numeric_list(payload, "alphas"),
        Alpha0s=_as_numeric_list(payload, "alpha0s"),
        Kappa0s=_as_numeric_list(payload, "kappa0s"),
        SYs=_as_numeric_list(payload, "sigma_ys"),
        Nrep=int(payload["n_rep"]),
        Njob=int(payload["n_jobs"]),
        Nrepfinal=int(payload["n_rep_final"]),
        fsel=float(payload["feature_selection_fraction"]),
    )

    proportions = pd.DataFrame(
        final_obj.ExpF(final_obj.Beta),
        index=mixture.columns,
        columns=mean.columns,
    )
    output_path = Path(payload["output_path"])
    output_path.parent.mkdir(parents=True, exist_ok=True)
    proportions.to_csv(output_path, sep="\t")

    if payload["extra_output_dir"] != "":
        extra_output_dir = Path(payload["extra_output_dir"])
        extra_output_dir.mkdir(parents=True, exist_ok=True)
        pd.DataFrame([best_set]).to_csv(extra_output_dir / "best_set.tsv", sep="\t", index=False)
        pd.DataFrame({"run": list(all_out.keys())}).to_csv(
            extra_output_dir / "all_runs.tsv",
            sep="\t",
            index=False,
        )
        pd.DataFrame({"best_e_step": [best_obj.E_step(best_obj.Nu, best_obj.Beta, best_obj.Omega)]}).to_csv(
            extra_output_dir / "best_score.tsv",
            sep="\t",
            index=False,
        )


if __name__ == "__main__":
    main()
