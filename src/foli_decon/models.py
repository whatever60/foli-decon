"""Data models used by the foli-decon package."""

from dataclasses import dataclass
from dataclasses import field

import pandas as pd


@dataclass
class DeconvolutionResult:
    """Container for normalized deconvolution outputs."""

    tool: str
    proportion: pd.DataFrame | None = None
    score: pd.DataFrame | None = None
    component: pd.DataFrame | None = None
    p_value: pd.DataFrame | None = None
    uncertainty: pd.DataFrame | None = None
    selected_features: pd.DataFrame | None = None
    metadata: dict[str, str | int | float | bool] = field(default_factory=dict)
    output_paths: dict[str, str] = field(default_factory=dict)


@dataclass
class ReferenceTrainingResult:
    """Container for trained reference artifacts."""

    tool: str
    signature: pd.DataFrame | None = None
    signature_path: str | None = None
    model_path: str | None = None
    model_dir: str | None = None
    rds_path: str | None = None
    h5ad_path: str | None = None
    metadata: dict[str, str | int | float | bool] = field(default_factory=dict)
    output_paths: dict[str, str] = field(default_factory=dict)


@dataclass
class FeatureSelectionResult:
    """Container for selected feature panels and ranking diagnostics."""

    tool: str
    selected_genes: pd.Index
    ranking: pd.DataFrame | None = None
    metadata: dict[str, str | int | float | bool] = field(default_factory=dict)
    output_paths: dict[str, str] = field(default_factory=dict)
