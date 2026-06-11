#!/usr/bin/env python
"""Download ABIS Shiny-app signature resources for local wrapper use."""

from pathlib import Path
from urllib.request import urlopen


ABIS_RAW_BASE_URL = "https://raw.githubusercontent.com/giannimonaco/ABIS/master/data"
OUTPUT_ROOT = Path(__file__).resolve().parents[1] / "resources" / "abis"
RESOURCES = {
    "rnaseq/sigmatrixRNAseq.txt": "sigmatrixRNAseq.txt",
    "microarray/sigmatrixMicro.txt": "sigmatrixMicro.txt",
    "microarray/target.txt": "target.txt",
}


def download_resource(relative_output: str, source_name: str) -> None:
    """Download one ABIS resource file into the repository resource folder."""

    output_path = OUTPUT_ROOT / relative_output
    output_path.parent.mkdir(parents=True, exist_ok=True)
    url = f"{ABIS_RAW_BASE_URL}/{source_name}"
    output_path.write_bytes(urlopen(url, timeout=60).read())
    print(f"{url} -> {output_path}")


def main() -> None:
    """Download all ABIS resources needed by the wrapper defaults."""

    for relative_output, source_name in RESOURCES.items():
        download_resource(relative_output=relative_output, source_name=source_name)


if __name__ == "__main__":
    main()
