#!/usr/bin/env python3
"""Run Protenix template/MSA sweeps for nanobody validation sets."""

from __future__ import annotations

import argparse
import csv
import json
import subprocess
import sys
import time
from pathlib import Path
from typing import Any

REPO_ROOT = Path(__file__).resolve().parents[2]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))


def _strip_templates(obj: Any) -> Any:
    if isinstance(obj, dict):
        return {
            key: _strip_templates(value)
            for key, value in obj.items()
            if key not in {"templatesPath", "template_path", "templatePath"}
        }
    if isinstance(obj, list):
        return [_strip_templates(value) for value in obj]
    return obj


def _prepare_mode_json(base_json: Path, output_json: Path, mode: str) -> str:
    with open(base_json, "r", encoding="utf-8") as handle:
        data = json.load(handle)
    note = ""
    if mode == "no_template":
        data = _strip_templates(data)
        note = "templates removed from input JSON"
    elif mode in {"top1", "top5"}:
        note = "top-k template filtering is best-effort; JSON/HHR filtering is not modified by this wrapper"
    output_json.parent.mkdir(parents=True, exist_ok=True)
    with open(output_json, "w", encoding="utf-8") as handle:
        json.dump(data, handle, indent=2)
    return note


def _read_summary(path: Path) -> dict[str, Any]:
    if not path.exists():
        return {}
    with open(path, "r", encoding="utf-8") as handle:
        return json.load(handle)


def _annotate_metrics_outputs(csv_path: Path, json_path: Path, mode: str, cutoff: str) -> None:
    if csv_path.exists():
        with open(csv_path, "r", newline="", encoding="utf-8") as handle:
            reader = csv.DictReader(handle)
            rows = list(reader)
            fieldnames = list(reader.fieldnames or [])
        for column in ("template_mode", "template_max_release_date"):
            if column not in fieldnames:
                fieldnames.append(column)
        for row in rows:
            row["template_mode"] = mode
            row["template_max_release_date"] = cutoff
        with open(csv_path, "w", newline="", encoding="utf-8") as handle:
            writer = csv.DictWriter(handle, fieldnames=fieldnames)
            writer.writeheader()
            writer.writerows(rows)

    if json_path.exists():
        with open(json_path, "r", encoding="utf-8") as handle:
            summary = json.load(handle)
        summary["template_mode"] = mode
        summary["template_max_release_date"] = cutoff
        with open(json_path, "w", encoding="utf-8") as handle:
            json.dump(summary, handle, indent=2)


def build_arg_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Run no-template/default-template/top-k template sweeps.")
    parser.add_argument("--manifest", required=True, type=Path, help="Validation/test canonical manifest.")
    parser.add_argument("--base-input-json", required=True, type=Path, help="Base Protenix inference JSON.")
    parser.add_argument("--checkpoint", required=True, type=Path, help="Released or fine-tuned checkpoint.")
    parser.add_argument("--out-root", required=True, type=Path, help="Template sweep output root.")
    parser.add_argument("--modes", default="no_template,default_template,top1,top5", help="Comma-separated template modes.")
    parser.add_argument("--seeds", default="101,102,103", help="Comma-separated inference seeds.")
    parser.add_argument("--model-name", default="protenix-v2", help="Model name.")
    parser.add_argument("--n-sample", type=int, default=5, help="Diffusion samples per seed.")
    parser.add_argument("--n-step", type=int, default=200, help="Diffusion steps.")
    parser.add_argument("--n-cycle", type=int, default=10, help="Recycling cycles.")
    parser.add_argument("--sample-mode", choices=["rank0", "all", "oracle"], default="rank0", help="Evaluation sample selection.")
    parser.add_argument("--template-max-release-date", default="2022-01-01", help="Template release cutoff recorded for leakage control.")
    parser.add_argument("--dry-run", action="store_true", help="Write configs and commands without launching inference/evaluation.")
    return parser


def main() -> None:
    args = build_arg_parser().parse_args()
    sweep_id = time.strftime("sweep_%Y%m%d_%H%M%S")
    sweep_dir = args.out_root / sweep_id
    config_dir = sweep_dir / "configs"
    pred_dir = sweep_dir / "predictions"
    sweep_dir.mkdir(parents=True, exist_ok=True)

    summary_rows: list[dict[str, Any]] = []
    commands: list[str] = []
    for mode in [part.strip() for part in args.modes.split(",") if part.strip()]:
        mode_json = config_dir / f"{mode}.json"
        note = _prepare_mode_json(args.base_input_json, mode_json, mode)
        use_template = "false" if mode == "no_template" else "true"
        mode_pred_dir = pred_dir / mode
        mode_metrics = sweep_dir / f"metrics_{mode}.csv"
        mode_summary = sweep_dir / f"summary_{mode}.json"

        infer_cmd = [
            sys.executable,
            "scripts/nanobody/run_inference_checkpoint.py",
            "--checkpoint",
            str(args.checkpoint),
            "--model-name",
            args.model_name,
            "--input-json",
            str(mode_json),
            "--out-dir",
            str(mode_pred_dir),
            "--seeds",
            args.seeds,
            "--n-sample",
            str(args.n_sample),
            "--n-step",
            str(args.n_step),
            "--n-cycle",
            str(args.n_cycle),
            "--use-template",
            use_template,
        ]
        eval_cmd = [
            sys.executable,
            "scripts/nanobody/evaluate_nanobody_predictions.py",
            "--manifest",
            str(args.manifest),
            "--pred-dir",
            str(mode_pred_dir),
            "--out-csv",
            str(mode_metrics),
            "--out-json",
            str(mode_summary),
            "--sample-mode",
            args.sample_mode,
        ]
        commands.extend([" ".join(infer_cmd), " ".join(eval_cmd)])
        if not args.dry_run:
            subprocess.run(infer_cmd, check=True)
            subprocess.run(eval_cmd, check=True)
            _annotate_metrics_outputs(
                mode_metrics,
                mode_summary,
                mode,
                args.template_max_release_date,
            )
        summary = _read_summary(mode_summary)
        summary_rows.append(
            {
                "sweep_id": sweep_id,
                "template_mode": mode,
                "template_max_release_date": args.template_max_release_date,
                "mean_dockq": summary.get("mean_dockq"),
                "dockq_success_rate_0_23": summary.get("dockq_success_rate_0_23"),
                "median_cdrh3_rmsd_ca": summary.get("median_cdrh3_rmsd_ca"),
                "mean_cdrh3_rmsd_ca": summary.get("mean_cdrh3_rmsd_ca"),
                "notes": note,
            }
        )

    with open(sweep_dir / "commands.sh", "w", encoding="utf-8") as handle:
        handle.write("\n".join(commands) + "\n")
    with open(sweep_dir / "summary.csv", "w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(summary_rows[0].keys()) if summary_rows else [])
        if summary_rows:
            writer.writeheader()
            writer.writerows(summary_rows)
    with open(sweep_dir / "summary.json", "w", encoding="utf-8") as handle:
        json.dump({"sweep_id": sweep_id, "modes": summary_rows}, handle, indent=2)
    print(sweep_dir)


if __name__ == "__main__":
    main()
