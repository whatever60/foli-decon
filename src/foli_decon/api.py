"""Unified high-level API for tool dispatch."""

from foli_decon.models import DeconvolutionResult, FeatureSelectionResult, ReferenceTrainingResult
from foli_decon.tools import (
    run_autogenes,
    run_abis,
    run_blade,
    run_cellcode,
    run_consensus_tme,
    run_bayesprism,
    run_bisque,
    run_blue,
    run_cdseq,
    run_cibersortx,
    run_deconrnaseq,
    run_dissect,
    run_dtangle,
    run_dwls,
    run_epic,
    run_estimate,
    run_instaprism,
    run_mead,
    run_mcp_counter,
    run_music,
    run_music2,
    run_psea,
    run_quantiseq,
    run_scdc,
    run_scaden,
    run_tape,
    run_xcell,
    run_xcell2,
    select_features_autogenes,
    select_features_markermap,
    select_features_scgenefit,
    train_xcell2_reference,
)


_UNSUPPORTED_BUT_REQUESTED_TOOLS: dict[str, str] = {
    "timer": "TIMER is blocked: local native support is unresolved, and web-only TIMER2/TIMER3 are not used by foli-decon.",
    "timer2": "TIMER2 is web-based and is not supported by the local foli-decon wrappers.",
    "timer3": "TIMER3 is web-based and is not supported by the local foli-decon wrappers.",
    "dsa": "DSA is blocked: no conda deconvolution package route was found; CRAN/PyPI name matches are unrelated packages.",
}

_TOOL_ALIASES: dict[str, str] = {
    "consensustme": "consensus_tme",
    "mcpcounter": "mcp_counter",
    "instaprime": "instaprism",
}


def supported_tools() -> tuple[str, ...]:
    """Return all supported deconvolution tool names."""

    return (
        "xcell",
        "quantiseq",
        "consensus_tme",
        "mcp_counter",
        "psea",
        "abis",
        "estimate",
        "autogenes",
        "blade",
        "blue",
        "cellcode",
        "cdseq",
        "xcell2",
        "epic",
        "dtangle",
        "bayesprism",
        "deconrnaseq",
        "cibersortx",
        "bisque",
        "dwls",
        "music",
        "music2",
        "mead",
        "scdc",
        "instaprism",
        "tape",
        "scaden",
        "dissect",
    )


def supported_reference_trainers() -> tuple[str, ...]:
    """Return all supported reference-training tool names."""

    return ("xcell2",)


def supported_feature_selectors() -> tuple[str, ...]:
    """Return all supported feature-selection tool names."""

    return ("autogenes", "markermap", "scgenefit")


def run_deconvolution(tool: str, **kwargs) -> DeconvolutionResult:
    """Dispatch to a specific deconvolution backend by name."""

    normalized_tool = tool.lower().replace("-", "_").replace(" ", "_")
    if normalized_tool in _TOOL_ALIASES:
        normalized_tool = _TOOL_ALIASES[normalized_tool]

    if normalized_tool == "xcell":
        return run_xcell(**kwargs)
    if normalized_tool == "quantiseq":
        return run_quantiseq(**kwargs)
    if normalized_tool in _UNSUPPORTED_BUT_REQUESTED_TOOLS:
        raise ValueError(_UNSUPPORTED_BUT_REQUESTED_TOOLS[normalized_tool])
    if normalized_tool == "consensus_tme":
        return run_consensus_tme(**kwargs)
    if normalized_tool == "mcp_counter":
        return run_mcp_counter(**kwargs)
    if normalized_tool == "psea":
        return run_psea(**kwargs)
    if normalized_tool == "abis":
        return run_abis(**kwargs)
    if normalized_tool == "estimate":
        return run_estimate(**kwargs)
    if normalized_tool == "autogenes":
        return run_autogenes(**kwargs)
    if normalized_tool == "blade":
        return run_blade(**kwargs)
    if normalized_tool == "blue":
        return run_blue(**kwargs)
    if normalized_tool == "cellcode":
        return run_cellcode(**kwargs)
    if normalized_tool == "cdseq":
        return run_cdseq(**kwargs)
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
    if normalized_tool == "music2":
        return run_music2(**kwargs)
    if normalized_tool == "mead":
        return run_mead(**kwargs)
    if normalized_tool == "scdc":
        return run_scdc(**kwargs)
    if normalized_tool == "instaprism":
        return run_instaprism(**kwargs)
    if normalized_tool == "tape":
        return run_tape(**kwargs)
    if normalized_tool == "scaden":
        return run_scaden(**kwargs)
    if normalized_tool == "dissect":
        return run_dissect(**kwargs)

    raise ValueError(f"Unsupported tool '{tool}'. Supported tools: {supported_tools()}")


def train_reference(tool: str, **kwargs) -> ReferenceTrainingResult:
    """Dispatch to a specific reference-training backend by name."""

    normalized_tool = tool.lower()

    if normalized_tool == "xcell2":
        return train_xcell2_reference(**kwargs)

    raise ValueError(f"Unsupported reference trainer '{tool}'. Supported trainers: {supported_reference_trainers()}")


def select_features(tool: str, **kwargs) -> FeatureSelectionResult:
    """Dispatch to a specific feature-selection backend by name."""

    normalized_tool = tool.lower().replace("-", "_").replace(" ", "_")
    if normalized_tool == "autogenes":
        return select_features_autogenes(**kwargs)
    if normalized_tool == "markermap":
        return select_features_markermap(**kwargs)
    if normalized_tool in {"scgenefit", "sc_gene_fit"}:
        return select_features_scgenefit(**kwargs)

    raise ValueError(f"Unsupported feature selector '{tool}'. Supported selectors: {supported_feature_selectors()}")
