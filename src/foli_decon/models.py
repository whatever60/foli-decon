"""Data models used by the foli-decon package."""

from dataclasses import dataclass

import pandas as pd


@dataclass
class DeconvolutionResult:
    """Container for normalized deconvolution outputs."""

    tool: str
    proportions: pd.DataFrame
    metadata: dict[str, str | int | float | bool]
