#!/usr/bin/env python3
"""Run Protenix inference from an arbitrary checkpoint file."""

from __future__ import annotations

import argparse
import os
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[2]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))


def _bool_string(value: str | bool) -> str:
    if isinstance(value, bool):
        return "true" if value else "false"
    return str(value).lower()


def build_arg_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Evaluate a released or fine-tuned checkpoint with runner/inference.py.")
    parser.add_argument("--checkpoint", required=True, type=Path, help="Checkpoint file to evaluate.")
    parser.add_argument("--model-name", default="protenix-v2", help="Model name; temp checkpoint will be named <model-name>.pt.")
    parser.add_argument("--input-json", required=True, type=Path, help="Protenix inference input JSON.")
    parser.add_argument("--out-dir", required=True, type=Path, help="Prediction output directory.")
    parser.add_argument("--seeds", default="101", help="Comma-separated inference seeds.")
    parser.add_argument("--n-sample", type=int, default=5, help="Diffusion samples per seed.")
    parser.add_argument("--n-step", type=int, default=200, help="Diffusion steps.")
    parser.add_argument("--n-cycle", type=int, default=10, help="Recycling cycles.")
    parser.add_argument("--use-msa", default="true", help="Whether to use MSA features.")
    parser.add_argument("--use-template", default="true", help="Whether to use templates.")
    parser.add_argument("--dtype", default="bf16", help="Inference dtype.")
    parser.add_argument("--triangle-attention", default="cuequivariance", help="Triangle attention kernel.")
    parser.add_argument("--triangle-multiplicative", default="cuequivariance", help="Triangle multiplicative kernel.")
    parser.add_argument("--copy-checkpoint", action="store_true", help="Copy instead of symlinking into the temporary checkpoint directory.")
    parser.add_argument("--extra-arg", action="append", default=[], help="Additional raw '--key value' argument pair for runner/inference.py.")
    return parser


def main() -> None:
    args = build_arg_parser().parse_args()
    checkpoint = args.checkpoint.expanduser().resolve()
    if not checkpoint.exists():
        raise FileNotFoundError(f"Checkpoint not found: {checkpoint}")
    args.out_dir.mkdir(parents=True, exist_ok=True)

    with tempfile.TemporaryDirectory(prefix="protenix_nanobody_ckpt_") as tmp:
        tmpdir = Path(tmp)
        linked = tmpdir / f"{args.model_name}.pt"
        if args.copy_checkpoint:
            shutil.copy2(checkpoint, linked)
        else:
            os.symlink(checkpoint, linked)

        cmd = [
            sys.executable,
            "runner/inference.py",
            "--model_name",
            args.model_name,
            "--load_checkpoint_dir",
            str(tmpdir),
            "--input_json_path",
            str(args.input_json),
            "--dump_dir",
            str(args.out_dir),
            "--seeds",
            args.seeds,
            "--sample_diffusion.N_sample",
            str(args.n_sample),
            "--sample_diffusion.N_step",
            str(args.n_step),
            "--model.N_cycle",
            str(args.n_cycle),
            "--use_msa",
            _bool_string(args.use_msa),
            "--use_template",
            _bool_string(args.use_template),
            "--dtype",
            args.dtype,
            "--triangle_attention",
            args.triangle_attention,
            "--triangle_multiplicative",
            args.triangle_multiplicative,
        ]
        for extra in args.extra_arg:
            cmd.extend(extra.split())
        raise SystemExit(subprocess.run(cmd, check=False).returncode)


if __name__ == "__main__":
    main()
