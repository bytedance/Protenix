#!/usr/bin/env python3
"""
Test script for RunPod Protenix API.

Supports all parameters from the updated runpod_handler.py
"""
import requests
import json
import os
import sys
import argparse
from pathlib import Path

try:
    from dotenv import load_dotenv
    load_dotenv()
except ImportError:
    pass  # dotenv not required


def parse_args():
    parser = argparse.ArgumentParser(description="Test RunPod Protenix API")
    parser.add_argument("--sync", action="store_true", default=True, 
                        help="Use synchronous endpoint (default)")
    parser.add_argument("--async", dest="sync", action="store_false",
                        help="Use asynchronous endpoint")
    parser.add_argument("--timeout", type=int, default=600,
                        help="Request timeout in seconds (default: 600)")
    parser.add_argument("--use-msa", type=bool, default=False,
                        help="Use MSA (default: False)")
    parser.add_argument("--seeds", type=int, nargs="+", default=[101],
                        help="Random seeds (default: [101])")
    parser.add_argument("--n-cycle", type=int, default=None,
                        help="Number of cycles (default: model-dependent)")
    parser.add_argument("--n-step", type=int, default=None,
                        help="Number of diffusion steps (default: model-dependent)")
    parser.add_argument("--n-sample", type=int, default=None,
                        help="Number of samples (default: 5)")
    parser.add_argument("--model", type=str, default=None,
                        help="Model name (default: protenix_base_default_v0.5.0)")
    parser.add_argument("--output-dir", type=str, default="./output",
                        help="Directory to save output files (default: ./output)")
    parser.add_argument("--sequence", type=str, default=None,
                        help="Protein sequence to predict (uses default test sequence if not provided)")
    return parser.parse_args()


def main():
    args = parse_args()
    
    # Configuration from environment
    API_KEY = os.environ.get("RUNPOD_API_KEY")
    ENDPOINT_ID = os.environ.get("RUNPOD_ENDPOINT_ID")
    
    if not API_KEY or API_KEY == "YOUR_API_KEY":
        print("ERROR: Set RUNPOD_API_KEY environment variable")
        sys.exit(1)
    if not ENDPOINT_ID or ENDPOINT_ID == "YOUR_ENDPOINT_ID":
        print("ERROR: Set RUNPOD_ENDPOINT_ID environment variable")
        sys.exit(1)
    
    # Choose endpoint
    endpoint_type = "runsync" if args.sync else "run"
    URL = f"https://api.runpod.ai/v2/{ENDPOINT_ID}/{endpoint_type}"
    
    headers = {
        "Content-Type": "application/json",
        "Authorization": f"Bearer {API_KEY}"
    }
    
    # Default test sequence if none provided
    default_sequence = "MAEVIRSSAFWRSFPIFEEFDSETLCELSGIASYRKWSAGTVIFQRGDQGDYMIVVVSGRIKLSLFTPQGRELMLRQHEAGALFGEMALLDGQPRSADATAVTAAEGYVIGKKDFLALITQRPKTAEAVIRFLCAQLRDTTDRLETIALYDLNARVARFFLATLRQIHGSEMPQSANLRLTLSQTDIASILGASRPKVNRAILSLEESGAIKRADGIICCNVGRLLSIADPEEDLEHHHHHHHH"
    sequence = args.sequence or default_sequence
    
    # Build payload
    payload = {
        "input": {
            "name": "test_request",
            "sequences": [
                {
                    "proteinChain": {
                        "sequence": sequence,
                        "count": 1
                    }
                }
            ],
            "use_msa": args.use_msa,
            "seeds": args.seeds,
        }
    }
    
    # Add optional parameters if provided
    if args.model:
        payload["input"]["model_name"] = args.model
    if args.n_cycle is not None:
        payload["input"]["n_cycle"] = args.n_cycle
    if args.n_step is not None:
        payload["input"]["n_step"] = args.n_step
    if args.n_sample is not None:
        payload["input"]["n_sample"] = args.n_sample
    
    print("=" * 60)
    print("RUNPOD PROTENIX API TEST")
    print("=" * 60)
    print(f"Endpoint: {URL}")
    print(f"Timeout: {args.timeout}s")
    print(f"Parameters:")
    print(f"  use_msa: {args.use_msa}")
    print(f"  seeds: {args.seeds}")
    print(f"  n_cycle: {args.n_cycle or 'default'}")
    print(f"  n_step: {args.n_step or 'default'}")
    print(f"  n_sample: {args.n_sample or 'default'}")
    print(f"  model: {args.model or 'default'}")
    print(f"  sequence length: {len(sequence)}")
    print("=" * 60)
    
    print(f"\nSending request...")
    try:
        response = requests.post(
            URL, 
            headers=headers, 
            json=payload,
            timeout=args.timeout
        )
        response.raise_for_status()
        
        result = response.json()
        print(f"Response Status Code: {response.status_code}")
        print(f"Response keys: {list(result.keys())}")
        
        # Check for errors in response
        if "error" in result:
            print("\n" + "=" * 60)
            print("ERROR IN RESPONSE")
            print("=" * 60)
            print(f"Error: {result.get('error')}")
            if "traceback" in result:
                print(f"\nTraceback:\n{result.get('traceback')}")
            return
        
        # Handle output
        if "output" in result:
            output_data = result["output"]
            
            # Check for error in output
            if "error" in output_data:
                print("\n" + "=" * 60)
                print("ERROR IN OUTPUT")
                print("=" * 60)
                print(f"Error: {output_data.get('error')}")
                if "traceback" in output_data:
                    print(f"\nTraceback:\n{output_data.get('traceback')}")
                return
            
            # Handle successful output with files
            if "output_files" in output_data:
                print("\n" + "=" * 60)
                print("SUCCESS! OUTPUT FILES RECEIVED")
                print("=" * 60)
                
                output_files = output_data["output_files"]
                print(f"Received {len(output_files)} file(s)")
                
                # Create output directory
                output_dir = Path(args.output_dir)
                output_dir.mkdir(parents=True, exist_ok=True)
                
                # Save files
                for fname, content in output_files.items():
                    # Create subdirectories if needed
                    fpath = output_dir / fname
                    fpath.parent.mkdir(parents=True, exist_ok=True)
                    
                    print(f"  Saving: {fname} ({len(content)} bytes)")
                    with open(fpath, "w") as f:
                        f.write(content)
                
                print(f"\nFiles saved to: {output_dir.absolute()}")
            else:
                print("\nOutput (no files):")
                print(json.dumps(output_data, indent=2))
        else:
            print("\nFull Response:")
            print(json.dumps(result, indent=2))
    
    except requests.exceptions.Timeout:
        print(f"\nERROR: Request timed out after {args.timeout}s")
        print("Consider using --async for long-running jobs")
    except requests.exceptions.RequestException as e:
        print(f"\nERROR: {e}")
        if hasattr(e, 'response') and e.response is not None:
            print(f"Response text: {e.response.text}")


if __name__ == "__main__":
    main()
