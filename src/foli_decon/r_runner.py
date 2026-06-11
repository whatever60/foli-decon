"""Execution helper for R-backed wrappers."""

import json
import os
import signal
import subprocess
import tempfile
from pathlib import Path

import pandas as pd

from foli_decon.io import read_sample_by_celltype


def get_r_script_path() -> Path:
    """Return the absolute path to the package R entrypoint script."""

    return Path(__file__).resolve().parent / "r_scripts" / "run_tools.R"


def run_r_tool(
    tool: str,
    payload: dict[str, object],
    rscript_command: list[str] | None = None,
    output_path: str | None = None,
    log_path: str | None = None,
) -> tuple[pd.DataFrame, dict[str, object]]:
    """Run an R method wrapper and return parsed output with enriched metadata."""

    with tempfile.TemporaryDirectory(prefix="foli_decon_r_") as temp_dir:
        temp_path = Path(temp_dir)
        payload_path = temp_path / "payload.json"
        resolved_output_path = temp_path / "output.tsv" if output_path is None else Path(output_path)
        resolved_output_path.parent.mkdir(parents=True, exist_ok=True)

        payload_copy = payload.copy()
        payload_copy["output_path"] = str(resolved_output_path)

        payload_path.write_text(json.dumps(payload_copy), encoding="utf-8")
        command = ["Rscript"] if rscript_command is None else rscript_command.copy()
        command.extend([str(get_r_script_path()), tool, str(payload_path)])
        timeout_raw = os.environ.get("FOLI_DECON_R_TIMEOUT_SECONDS", "")
        timeout_seconds = None if timeout_raw == "" else float(timeout_raw)
        process = subprocess.Popen(
            command,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
            start_new_session=True,
        )
        try:
            stdout, stderr = process.communicate(timeout=timeout_seconds)
        except subprocess.TimeoutExpired as error:
            os.killpg(process.pid, signal.SIGTERM)
            try:
                stdout, stderr = process.communicate(timeout=10)
            except subprocess.TimeoutExpired:
                os.killpg(process.pid, signal.SIGKILL)
                stdout, stderr = process.communicate()
            error.output = stdout
            error.stderr = stderr
            raise
        if log_path is not None:
            resolved_log_path = Path(log_path)
            resolved_log_path.parent.mkdir(parents=True, exist_ok=True)
            resolved_log_path.write_text(
                f"stdout:\n{stdout}\n\nstderr:\n{stderr}",
                encoding="utf-8",
            )
        if process.returncode != 0:
            raise RuntimeError(
                f"R tool '{tool}' failed with exit code {process.returncode}.\n"
                f"stdout:\n{stdout}\n"
                f"stderr:\n{stderr}"
            )

        result = read_sample_by_celltype(str(resolved_output_path))
    return result, payload_copy
