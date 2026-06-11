"""Execution helpers for Python/CLI-backed deconvolution wrappers."""

import os
import subprocess
from pathlib import Path


def run_command(command: list[str], cwd: str | None = None, log_path: str | None = None) -> subprocess.CompletedProcess:
    """Run a subprocess command and raise with captured logs on failure."""

    timeout_raw = os.environ["FOLI_DECON_PY_TIMEOUT_SECONDS"] if "FOLI_DECON_PY_TIMEOUT_SECONDS" in os.environ else ""
    timeout_seconds = None if timeout_raw == "" else float(timeout_raw)
    process = subprocess.run(
        command,
        cwd=cwd,
        capture_output=True,
        text=True,
        timeout=timeout_seconds,
        check=False,
    )
    if log_path is not None:
        resolved_log_path = Path(log_path)
        resolved_log_path.parent.mkdir(parents=True, exist_ok=True)
        resolved_log_path.write_text(
            f"command:\n{' '.join(command)}\n\nstdout:\n{process.stdout}\n\nstderr:\n{process.stderr}",
            encoding="utf-8",
        )
    if process.returncode != 0:
        raise RuntimeError(
            f"Command failed with exit code {process.returncode}: {' '.join(command)}\n"
            f"stdout:\n{process.stdout}\n"
            f"stderr:\n{process.stderr}"
        )
    return process
