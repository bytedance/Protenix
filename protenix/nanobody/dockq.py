"""Optional DockQ integration for nanobody-antigen complex evaluation."""

from __future__ import annotations

import re
import shutil
import subprocess
from pathlib import Path
from typing import Any


def find_dockq_executable() -> str | None:
    for name in ("DockQ", "DockQ.py"):
        exe = shutil.which(name)
        if exe:
            return exe
    return None


def parse_dockq_output(output: str) -> dict[str, float | None]:
    patterns = {
        "dockq": r"DockQ\s*[:=]?\s*([0-9]*\.?[0-9]+)",
        "irmsd": r"iRMSD\s*[:=]?\s*([0-9]*\.?[0-9]+)",
        "lrmsd": r"LRMSD\s*[:=]?\s*([0-9]*\.?[0-9]+)",
        "fnat": r"Fnat\s*[:=]?\s*([0-9]*\.?[0-9]+)",
    }
    parsed: dict[str, float | None] = {}
    for key, pattern in patterns.items():
        match = re.search(pattern, output, flags=re.IGNORECASE)
        parsed[key] = float(match.group(1)) if match else None
    return parsed


def compute_dockq(
    pred_path: str | Path,
    native_path: str | Path,
    model_chains: str | None = None,
    native_chains: str | None = None,
    timeout: int = 300,
) -> dict[str, Any]:
    """Run the optional DockQ CLI if available."""

    executable = find_dockq_executable()
    if executable is None:
        return {
            "dockq_available": False,
            "dockq": None,
            "dockq_success_0_23": None,
            "irmsd": None,
            "lrmsd": None,
            "fnat": None,
            "metric_status": "dockq_unavailable",
            "notes": "DockQ executable not found; install the optional DockQ package to enable complex metrics.",
        }

    cmd = [executable, str(pred_path), str(native_path)]
    if model_chains:
        cmd.extend(["--model_chain", model_chains])
    if native_chains:
        cmd.extend(["--native_chain", native_chains])

    try:
        proc = subprocess.run(
            cmd,
            text=True,
            capture_output=True,
            check=False,
            timeout=timeout,
        )
    except Exception as exc:  # pragma: no cover - depends on optional CLI
        return {
            "dockq_available": True,
            "dockq": None,
            "dockq_success_0_23": None,
            "irmsd": None,
            "lrmsd": None,
            "fnat": None,
            "metric_status": "dockq_failed",
            "notes": str(exc),
        }

    output = "\n".join([proc.stdout, proc.stderr])
    parsed = parse_dockq_output(output)
    dockq = parsed.get("dockq")
    return {
        "dockq_available": dockq is not None,
        "dockq": dockq,
        "dockq_success_0_23": (dockq > 0.23) if dockq is not None else None,
        "irmsd": parsed.get("irmsd"),
        "lrmsd": parsed.get("lrmsd"),
        "fnat": parsed.get("fnat"),
        "metric_status": "ok" if dockq is not None else "dockq_failed",
        "notes": "" if dockq is not None else output.strip(),
    }
