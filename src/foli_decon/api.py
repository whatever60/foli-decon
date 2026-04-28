"""Unified high-level API for tool dispatch."""

from foli_decon.models import DeconvolutionResult
from foli_decon.tools import (
    run_bayesprism,
    run_bisque,
    run_cibersortx,
    run_deconrnaseq,
    run_dtangle,
    run_dwls,
    run_epic,
    run_music,
    run_xcell,
    run_xcell2,
)


def supported_tools() -> tuple[str, ...]:
    """Return all supported deconvolution tool names."""

    return (
        "xcell",
        "xcell2",
        "epic",
        "dtangle",
        "bayesprism",
        "deconrnaseq",
        "cibersortx",
        "bisque",
        "dwls",
        "music",
    )


def run_deconvolution(tool: str, **kwargs) -> DeconvolutionResult:
    """Dispatch to a specific deconvolution backend by name."""

    normalized_tool = tool.lower()

    if normalized_tool == "xcell":
        return run_xcell(**kwargs)
    if normalized_tool == "xcell2":
        return run_xcell2(**kwargs)
    if normalized_tool == "epic":
        return run_epic(**kwargs)
    if normalized_tool == "dtangle":
        return run_dtangle(**kwargs)
    if normalized_tool == "bayesprism":
        return run_bayesprism(**kwargs)
    if normalized_tool == "deconrnaseq":
        return run_deconrnaseq(**kwargs)
    if normalized_tool == "cibersortx":
        return run_cibersortx(**kwargs)
    if normalized_tool == "bisque":
        return run_bisque(**kwargs)
    if normalized_tool == "dwls":
        return run_dwls(**kwargs)
    if normalized_tool == "music":
        return run_music(**kwargs)

    raise ValueError(f"Unsupported tool '{tool}'. Supported tools: {supported_tools()}")
