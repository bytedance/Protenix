"""
RunPod Serverless Handler for Protenix Inference.

Uses the proven infer_predict() pipeline approach that works in the non-serverless deployment.
"""
import logging
import os
import shutil
import tempfile
import json
import runpod
import torch
import traceback
from argparse import Namespace
from typing import Any, Dict
from ml_collections.config_dict import ConfigDict

# Import Protenix modules
from configs.configs_base import configs as configs_base
from configs.configs_data import data_configs
from configs.configs_inference import inference_configs
from configs.configs_model_type import model_configs
from protenix.config import parse_configs
from runner.inference import (
    InferenceRunner,
    infer_predict,
    update_gpu_compatible_configs,
    download_infercence_cache,
)

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)

# PyTorch 2.6+ compatibility fix for ESM model loading
# The fair-esm repository is archived and can't be updated to support newer PyTorch versions.
# ESM model files contain argparse.Namespace which isn't allowed by default in secure unpickling.
torch.serialization.add_safe_globals([Namespace])

# Global configs (initialized once on cold start)
base_configs = None


def init_handler():
    """
    Initialize base configurations and download checkpoints once during container start/cold start.
    """
    global base_configs
    
    logger.info("Initializing handler...")
    
    # Setup base configs - merge configs similar to runner/inference.py
    configs = {**configs_base, **{"data": data_configs}, **inference_configs}
    
    # Set a placeholder for input_json_path (will be overridden per request)
    configs["input_json_path"] = "/dev/null"
    
    merged_configs = ConfigDict(configs)
    
    # Set default values usually coming from CLI
    merged_configs.model_name = os.environ.get("MODEL_NAME", "protenix_base_default_v0.5.0")
    merged_configs.dump_dir = os.environ.get("DUMP_DIR", "/tmp/output")
    merged_configs.load_checkpoint_dir = os.environ.get("CHECKPOINT_DIR", "/app/checkpoints")
    merged_configs.need_atom_confidence = False
    merged_configs.use_msa = True  # Default true, can be overridden in payload
    merged_configs.num_workers = 0
    merged_configs.dtype = "bf16"
    merged_configs.seeds = [101]  # Default seed
    merged_configs.sample_diffusion = ConfigDict({
        "step_scale_eta": 1.5,
        "N_sample": 1,
        "N_step": 50
    })
    merged_configs.enable_tf32 = True
    merged_configs.enable_efficient_fusion = True
    merged_configs.enable_diffusion_shared_vars_cache = True
    merged_configs.infer_setting = ConfigDict({
        "sample_diffusion_chunk_size": 5,
        "chunk_size": 256
    })
    merged_configs.triangle_multiplicative = "torch"
    merged_configs.load_strict = True
    merged_configs.skip_amp = ConfigDict({
        "confidence_head": False,
        "sample_diffusion": False
    })
    merged_configs.sorted_by_ranking_score = True
    merged_configs.deterministic = False
    
    # Update model specific configs
    model_name = merged_configs.model_name
    if model_name in model_configs:
        merged_configs.update(ConfigDict(model_configs[model_name]))
    
    merged_configs = update_gpu_compatible_configs(merged_configs)
    
    # Download checkpoints/cache if needed (do this once on cold start)
    logger.info("Downloading inference cache if needed...")
    download_infercence_cache(merged_configs)
    
    base_configs = merged_configs
    logger.info("Handler initialized successfully.")


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
            except UnicodeDecodeError:
                # Skip binary files
                logger.warning(f"Skipping binary file: {filepath}")
            except Exception as e:
                logger.error(f"Error reading file {filepath}: {e}")
    
    logger.info(f"Collected {len(output_files)} output files")
    return output_files


def handler(event):
    """
    RunPod Handler function.
    
    Uses the proven infer_predict() pipeline approach.
    """
    temp_dir = None
    
    try:
        job_input = event["input"]
        job_id = event.get("id", "unknown")
        
        # Ensure name exists
        if "name" not in job_input:
            job_input["name"] = f"job_{job_id}"
        
        sample_name = job_input["name"]
        logger.info(f"Processing job {job_id} for sample: {sample_name}")
        
        # Create temporary directory for this job
        temp_dir = tempfile.mkdtemp(prefix=f"protenix_{job_id}_")
        logger.info(f"Created temp directory: {temp_dir}")
        
        # Save input to JSON file
        json_path = save_input_to_json(job_input, temp_dir)
        
        # Create job-specific output directory
        output_dir = os.path.join(temp_dir, "output")
        os.makedirs(output_dir, exist_ok=True)
        
        # Create job-specific configs (copy from base)
        configs = ConfigDict(base_configs.to_dict())
        configs.input_json_path = json_path
        configs.dump_dir = output_dir
        
        # Override configs from request payload if provided
        if "use_msa" in job_input:
            configs.use_msa = job_input["use_msa"]
        if "seeds" in job_input:
            configs.seeds = job_input["seeds"] if isinstance(job_input["seeds"], list) else [job_input["seeds"]]
        if "n_sample" in job_input:
            configs.sample_diffusion.N_sample = job_input["n_sample"]
        if "n_step" in job_input:
            configs.sample_diffusion.N_step = job_input["n_step"]
        if "n_cycle" in job_input:
            configs.model.N_cycle = job_input["n_cycle"]
        
        logger.info(f"Running inference with: use_msa={configs.use_msa}, seeds={configs.seeds}")
        
        # Initialize runner and run inference using the proven pipeline
        runner = InferenceRunner(configs)
        infer_predict(runner, configs)
        
        # Clear CUDA cache
        if torch.cuda.is_available():
            torch.cuda.empty_cache()
        
        # Collect output files
        output_files = collect_output_files(output_dir)
        
        if not output_files:
            logger.warning("No output files found after inference")
            return {"error": "Inference completed but no output files were generated"}
        
        logger.info(f"Job {job_id} completed successfully with {len(output_files)} files")
        
        return {"output_files": output_files}
    
    except Exception as e:
        logger.error(f"Error processing job: {e}")
        logger.error(traceback.format_exc())
        
        # Clear CUDA cache on error
        if torch.cuda.is_available():
            torch.cuda.empty_cache()
        
        return {"error": str(e), "traceback": traceback.format_exc()}
    
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
