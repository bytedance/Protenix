"""
RunPod Serverless Handler for Protenix Inference.

Based on the working dropxcell-proteinx/api_server/inference_engine.py pattern.
"""
import logging
import os
import shutil
import tempfile
import json
import runpod
import torch
import traceback
import platform
import sys
from argparse import Namespace
from typing import Any, Dict, Optional

# Configure logging first
LOG_FORMAT = "%(asctime)s [%(filename)s:%(lineno)d] %(levelname)s %(name)s: %(message)s"
logging.basicConfig(
    level=logging.INFO,
    format=LOG_FORMAT,
    datefmt="%Y-%m-%d %H:%M:%S",
)
logger = logging.getLogger(__name__)

# Import Protenix modules
from ml_collections.config_dict import ConfigDict
from configs.configs_base import configs as configs_base
from configs.configs_data import data_configs
from configs.configs_inference import inference_configs
from configs.configs_model_type import model_configs
from protenix.config import parse_configs
from runner.inference import (
    InferenceRunner,
    infer_predict,
    download_infercence_cache,
)

# PyTorch 2.6+ compatibility fix for ESM model loading
# The fair-esm repository is archived and can't be updated to support newer PyTorch versions.
# ESM model files contain argparse.Namespace which isn't allowed by default in secure unpickling.
torch.serialization.add_safe_globals([Namespace])

# Global state (initialized once on cold start)
cache_downloaded = False


def log_separator(title: str = ""):
    """Log a visual separator for clarity."""
    logger.info("=" * 60)
    if title:
        logger.info(title)
        logger.info("=" * 60)


def log_runtime_info():
    """Log runtime environment information for debugging."""
    logger.info("RUNTIME ENVIRONMENT:")
    logger.info(f"  Python: {platform.python_version()}")
    logger.info(f"  Platform: {platform.platform()}")
    logger.info(f"  PyTorch: {torch.__version__}")
    logger.info(f"  CUDA available: {torch.cuda.is_available()}")
    if torch.cuda.is_available():
        logger.info(f"  CUDA device count: {torch.cuda.device_count()}")
        logger.info(f"  CUDA device name: {torch.cuda.get_device_name(0)}")


def init_handler():
    """
    Initialize handler - download checkpoints/cache once on cold start.
    
    This is called once when the container starts.
    """
    global cache_downloaded
    
    log_separator("INIT_HANDLER STARTING")
    log_runtime_info()
    
    # Create a minimal config just to download the cache
    logger.info("Creating minimal config for cache download...")
    
    local_inference_configs = dict(inference_configs)
    local_inference_configs["dump_dir"] = "/tmp/init_output"
    local_inference_configs["input_json_path"] = "/dev/null"
    local_inference_configs["model_name"] = os.environ.get("MODEL_NAME", "protenix_base_default_v0.5.0")
    
    configs = {**configs_base, **{"data": data_configs}, **local_inference_configs}
    configs = parse_configs(
        configs=configs,
        fill_required_with_null=True,
    )
    
    model_name = configs.model_name
    logger.info(f"Model name: {model_name}")
    
    # Update with model-specific configs
    if model_name in model_configs:
        model_specific_configs = ConfigDict(model_configs[model_name])
        configs.update(model_specific_configs)
        logger.info(f"Applied model-specific configs for: {model_name}")
    
    # Download checkpoints/cache if needed (do this once on cold start)
    logger.info("Downloading inference cache if needed...")
    download_infercence_cache(configs)
    cache_downloaded = True
    
    log_separator("INIT_HANDLER COMPLETED")


def save_input_to_json(input_data: Dict[str, Any], output_dir: str) -> str:
    """
    Save input payload to a JSON file in the format expected by Protenix.
    
    Args:
        input_data: The input dictionary from the RunPod request
        output_dir: Directory to save the JSON file
        
    Returns:
        Path to the saved JSON file
    """
    # The input_data should already be in Protenix format (with "name", "sequences", etc.)
    # Wrap in a list as Protenix expects a list of samples
    json_data = [input_data]
    
    json_path = os.path.join(output_dir, "input.json")
    with open(json_path, "w") as f:
        json.dump(json_data, f, indent=2)
    
    logger.info(f"Saved input JSON to: {json_path}")
    return json_path


def collect_output_files(output_dir: str) -> Dict[str, str]:
    """
    Collect all output files from the inference output directory.
    
    Args:
        output_dir: The output directory from inference
        
    Returns:
        Dictionary mapping filename to file contents
    """
    output_files = {}
    
    # Walk through all subdirectories to find output files
    for root, dirs, files in os.walk(output_dir):
        for filename in files:
            filepath = os.path.join(root, filename)
            try:
                # Try to read as text (CIF, JSON files are text)
                with open(filepath, "r") as f:
                    content = f.read()
                # Use relative path from output_dir as key
                rel_path = os.path.relpath(filepath, output_dir)
                output_files[rel_path] = content
                logger.info(f"  Collected: {rel_path} ({len(content)} bytes)")
            except UnicodeDecodeError:
                # Skip binary files
                logger.warning(f"  Skipping binary file: {filepath}")
            except Exception as e:
                logger.error(f"  Error reading file {filepath}: {e}")
    
    logger.info(f"Collected {len(output_files)} output files total")
    return output_files


def handler(event):
    """
    RunPod Handler function.
    
    Uses the proven pattern from dropxcell-proteinx inference_engine.py
    """
    temp_dir = None
    
    try:
        log_separator("HANDLER STARTED")
        
        job_input = event.get("input", {})
        job_id = event.get("id", "unknown")
        
        logger.info(f"Job ID: {job_id}")
        logger.info(f"Event keys: {list(event.keys())}")
        logger.info(f"Input keys: {list(job_input.keys())}")
        
        # Ensure name exists
        if "name" not in job_input:
            job_input["name"] = f"job_{job_id}"
        
        sample_name = job_input["name"]
        logger.info(f"Processing sample: {sample_name}")
        
        # Create temporary directory for this job
        temp_dir = tempfile.mkdtemp(prefix=f"protenix_{job_id}_")
        logger.info(f"Created temp directory: {temp_dir}")
        
        # Save input to JSON file
        json_path = save_input_to_json(job_input, temp_dir)
        
        # Create job-specific output directory
        output_dir = os.path.join(temp_dir, "output")
        os.makedirs(output_dir, exist_ok=True)
        logger.info(f"Output directory: {output_dir}")
        
        # Extract parameters from request (with defaults)
        model_name = job_input.get("model_name", os.environ.get("MODEL_NAME", "protenix_base_default_v0.5.0"))
        use_msa = job_input.get("use_msa", True)
        seeds = job_input.get("seeds", [101])
        n_cycle = job_input.get("n_cycle", None)
        n_step = job_input.get("n_step", None)
        n_sample = job_input.get("n_sample", None)
        
        # Ensure seeds is a list
        if not isinstance(seeds, list):
            seeds = [seeds]
        
        logger.info(f"Parameters:")
        logger.info(f"  model_name: {model_name}")
        logger.info(f"  use_msa: {use_msa}")
        logger.info(f"  seeds: {seeds}")
        logger.info(f"  n_cycle: {n_cycle}")
        logger.info(f"  n_step: {n_step}")
        logger.info(f"  n_sample: {n_sample}")
        
        # Set defaults based on model type (matching official Protenix implementation)
        if "mini" in model_name or "tiny" in model_name:
            default_n_cycle = 4
            default_n_step = 5
        else:
            default_n_cycle = 10
            default_n_step = 200
        
        # Handle None, 0, or any falsy value by using defaults
        n_cycle = n_cycle if (n_cycle is not None and n_cycle > 0) else default_n_cycle
        n_step = n_step if (n_step is not None and n_step > 0) else default_n_step
        n_sample = n_sample if (n_sample is not None and n_sample > 0) else 5
        
        logger.info(f"Final parameters after defaults:")
        logger.info(f"  n_cycle: {n_cycle}")
        logger.info(f"  n_step: {n_step}")
        logger.info(f"  n_sample: {n_sample}")
        
        # ============================================================
        # CONFIGURE INFERENCE - Following inference_engine.py pattern
        # ============================================================
        log_separator("CONFIGURING INFERENCE")
        
        # Use copy to avoid global mutation
        local_inference_configs = dict(inference_configs)
        local_inference_configs["dump_dir"] = output_dir
        local_inference_configs["input_json_path"] = json_path
        local_inference_configs["model_name"] = model_name
        
        logger.info("Merging configs...")
        configs = {**configs_base, **{"data": data_configs}, **local_inference_configs}
        
        logger.info("Calling parse_configs()...")
        configs = parse_configs(
            configs=configs,
            fill_required_with_null=True,
        )
        logger.info("parse_configs() completed successfully")
        
        # Update with model-specific configs
        if model_name in model_configs:
            model_specific_configs = ConfigDict(model_configs[model_name])
            configs.update(model_specific_configs)
            logger.info(f"Applied model-specific configs for: {model_name}")
        else:
            logger.warning(f"No model-specific configs found for: {model_name}")
        
        # Set user-provided parameters - PLAIN PYTHON TYPES, NOT ListValue
        logger.info("Setting user parameters (plain Python types)...")
        configs.seeds = seeds          # plain list
        configs.use_msa = use_msa      # plain bool
        configs.model.N_cycle = n_cycle
        configs.sample_diffusion.N_step = n_step
        configs.sample_diffusion.N_sample = n_sample
        
        logger.info(f"Config seeds type: {type(configs.seeds)}")
        logger.info(f"Config use_msa type: {type(configs.use_msa)}")
        
        # Download cache if needed (should be cached from init_handler)
        if not cache_downloaded:
            logger.info("Cache not downloaded yet, downloading now...")
            download_infercence_cache(configs)
        else:
            logger.info("Cache already downloaded during init")
        
        # ============================================================
        # RUN INFERENCE
        # ============================================================
        log_separator("INITIALIZING INFERENCE RUNNER")
        
        logger.info("Creating InferenceRunner...")
        runner = InferenceRunner(configs)
        logger.info("InferenceRunner created successfully!")
        
        logger.info("Starting infer_predict()...")
        infer_predict(runner, configs)
        logger.info("infer_predict() completed!")
        
        # Clear CUDA cache
        if torch.cuda.is_available():
            torch.cuda.empty_cache()
            logger.info("Cleared CUDA cache")
        
        # ============================================================
        # COLLECT RESULTS
        # ============================================================
        log_separator("COLLECTING RESULTS")
        
        output_files = collect_output_files(output_dir)
        
        if not output_files:
            logger.warning("No output files found after inference")
            return {"error": "Inference completed but no output files were generated"}
        
        log_separator(f"JOB {job_id} COMPLETED SUCCESSFULLY")
        logger.info(f"Returning {len(output_files)} files")
        
        return {"output_files": output_files}
    
    except Exception as e:
        logger.error(f"Error processing job: {e}")
        logger.error(traceback.format_exc())
        
        # Clear CUDA cache on error
        if torch.cuda.is_available():
            torch.cuda.empty_cache()
        
        return {
            "error": str(e),
            "traceback": traceback.format_exc(),
        }
    
    finally:
        # Cleanup temp directory
        if temp_dir and os.path.exists(temp_dir):
            try:
                shutil.rmtree(temp_dir)
                logger.info(f"Cleaned up temp directory: {temp_dir}")
            except Exception as e:
                logger.error(f"Error cleaning up temp directory: {e}")


if __name__ == "__main__":
    init_handler()
    runpod.serverless.start({"handler": handler})
