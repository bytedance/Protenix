# Copyright 2024 ByteDance and/or its affiliates.
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#      http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

import json
import os
import concurrent.futures
import multiprocessing
from functools import partial
from os.path import join as opjoin
from typing import Callable, Dict, Tuple, List, Set, Any, Mapping, Union, Optional
from tqdm import tqdm
import time
import fcntl  # For file locking

from utils import (
    convert_to_shared_dict,  # To create new shared dictionaries
    SharedDict,              # To handle type annotation
    release_shared_dict,     # To manually release dictionaries
    get_shared_dict_ids      # To list available dictionaries
)

# Type alias for dictionary-like objects (regular dict or Manager.dict)
DictLike = Union[Dict[str, Any], Mapping[str, Any], SharedDict]


def load_mapping_data(seq_to_pdb_id_path: str, seq_to_pdb_index_path: str, use_shared_memory: bool = False) -> Tuple[Dict[str, Any], DictLike, DictLike]:
    """
    Load mapping data from JSON files.
    
    Args:
        seq_to_pdb_id_path: Path to the seq_to_pdb_id_entity_id.json file
        seq_to_pdb_index_path: Path to the seq_to_pdb_index.json file
        use_shared_memory: Whether to use shared memory for dictionaries
        
    Returns:
        Tuple containing (seq_to_pdbid, first_pdbid_to_seq, seq_to_pdb_index) dictionaries
    """
    # Load sequence to PDB ID mapping
    with open(seq_to_pdb_id_path, "r") as f:
        seq_to_pdbid: Dict[str, Any] = json.load(f)
    
    # Create reverse mapping for easy lookup
    first_pdbid_to_seq_data = {"_".join(v[0]): k for k, v in seq_to_pdbid.items()}
    
    # Load sequence to PDB index mapping
    with open(seq_to_pdb_index_path, "r") as f:
        seq_to_pdb_index_data = json.load(f)
    
    # If using shared memory, convert the dictionaries to shared objects
    if use_shared_memory:
        # Create shared dictionaries
        first_pdbid_to_seq = convert_to_shared_dict(first_pdbid_to_seq_data)
        seq_to_pdb_index = convert_to_shared_dict(seq_to_pdb_index_data)
        
        print(f"Created shared memory dictionaries: {len(first_pdbid_to_seq)} PDB IDs, {len(seq_to_pdb_index)} index mappings")
    else:
        first_pdbid_to_seq = first_pdbid_to_seq_data
        seq_to_pdb_index = seq_to_pdb_index_data
        
    return seq_to_pdbid, first_pdbid_to_seq, seq_to_pdb_index


def rematch(pdb_line: str, first_pdbid_to_seq: DictLike, seq_to_pdb_index: DictLike) -> Tuple[str, str]:
    """
    Match a PDB line to its corresponding sequence and index.
    
    Args:
        pdb_line: PDB header line
        first_pdbid_to_seq: Dictionary mapping PDB IDs to sequences
        seq_to_pdb_index: Dictionary mapping sequences to PDB indices
        
    Returns:
        Tuple of (pdb_index, origin_query_seq)
    """
    pdb_id = pdb_line[1:-1]
    origin_query_seq = first_pdbid_to_seq[pdb_id]
    pdb_index = seq_to_pdb_index[origin_query_seq]
    return pdb_index, origin_query_seq


def write_log(
    msg: str,
    fname: str,
    log_root: str,
) -> None:
    """
    Write a log message to a file with proper file locking to handle concurrency.
    
    Args:
        msg: Message to log
        fname: File name associated with the log
        log_root: Root directory for log files
    """
    basename = fname.split(".")[0]
    log_path = opjoin(log_root, f"{basename}-{msg}")
    
    # Create a directory for lock files
    lock_dir = opjoin(log_root, "locks")
    os.makedirs(lock_dir, exist_ok=True)
    
    # Use a separate lock file for each log file
    lock_path = opjoin(lock_dir, f"{basename}-{msg}.lock")
    
    try:
        # Open (or create) the lock file
        with open(lock_path, 'w') as lock_file:
            # Acquire an exclusive lock (blocking)
            fcntl.flock(lock_file, fcntl.LOCK_EX)
            
            # Now safely create the log file
            with open(log_path, "w") as f:
                f.write(msg)
            
            # The lock is automatically released when the file is closed
    except Exception as e:
        # If something goes wrong, log it but don't crash
        print(f"Warning: Failed to write log for {fname}: {e}")


def process_one_file(
    fname: str, 
    msa_root: str, 
    save_root: str, 
    logger: Callable,
    first_pdbid_to_seq: DictLike,
    seq_to_pdb_index: DictLike
) -> None:
    """
    Process a single MSA file.
    
    Args:
        fname: Filename of the MSA file to process
        msa_root: Root directory containing MSA files
        save_root: Root directory to save processed files
        logger: Function to log events
        first_pdbid_to_seq: Dictionary mapping PDB IDs to sequences
        seq_to_pdb_index: Dictionary mapping sequences to PDB indices
    """
    with open(file_path := opjoin(msa_root, fname), "r") as f:
        for i, line in enumerate(f):
            if i == 0:
                pdb_line = line
            if i == 1:
                if len(line) == 1:
                    logger("empty_query_seq", fname)
                    return
                query_line = line
                break

    save_fname, origin_query_seq = rematch(pdb_line, first_pdbid_to_seq, seq_to_pdb_index)

    os.makedirs(sub_dir_path := opjoin(save_root, f"{save_fname}"), exist_ok=True)
    uniref100_lines = [">query\n", f"{origin_query_seq}\n"]
    other_lines = [">query\n", f"{origin_query_seq}\n"]

    with open(file_path, "r") as f:
        lines = f.readlines()

    for i, line in enumerate(lines):
        if i < 2:
            continue
        if i % 2 == 0:
            # header
            if not line.startswith(">"):
                logger(f"bad_header_{i}", fname)
                return
            seq = lines[i + 1]

            if line.startswith(">UniRef100"):
                uniref100_lines.extend([line, seq])
            else:
                other_lines.extend([line, seq])

    assert len(other_lines) + len(uniref100_lines) - 2 == len(lines)

    other_lines = other_lines[0:2] + other_lines[4:]
    for i, line in enumerate(other_lines):
        if i > 0 and i % 2 == 0:
            assert "\t" in line
    with open(opjoin(sub_dir_path, "uniref100_hits.a3m"), "w") as f:
        for line in uniref100_lines:
            f.write(line)
    with open(opjoin(sub_dir_path, "mmseqs_other_hits.a3m"), "w") as f:
        for line in other_lines:
            f.write(line)


def process_file_batch(
    file_batch: List[str], 
    msa_root: str, 
    save_root: str, 
    log_root: str,
    first_pdbid_to_seq: DictLike,
    seq_to_pdb_index: DictLike
) -> Set[str]:
    """
    Process a batch of MSA files.
    
    Args:
        file_batch: List of filenames to process in this batch
        msa_root: Root directory containing MSA files
        save_root: Root directory to save processed files
        log_root: Root directory for log files
        first_pdbid_to_seq: Dictionary mapping PDB IDs to sequences
        seq_to_pdb_index: Dictionary mapping sequences to PDB indices
        
    Returns:
        Set of files that were processed successfully
    """
    # Create a logger for this batch
    batch_logger = partial(write_log, log_root=log_root)
    
    # Track completed files
    completed_files = set()
    
    # Process each file in the batch
    for fname in file_batch:
        try:
            process_one_file(
                fname=fname,
                msa_root=msa_root,
                save_root=save_root,
                logger=batch_logger,
                first_pdbid_to_seq=first_pdbid_to_seq,
                seq_to_pdb_index=seq_to_pdb_index
            )
            completed_files.add(fname)
        except Exception as e:
            # Log any exceptions but continue processing the batch
            basename = fname.split(".")[0]
            error_path = opjoin(log_root, f"{basename}-exception")
            lock_dir = opjoin(log_root, "locks")
            lock_path = opjoin(lock_dir, f"{basename}-exception.lock")
            
            try:
                # Ensure lock directory exists
                os.makedirs(lock_dir, exist_ok=True)
                
                # Use file locking for the error log too
                with open(lock_path, 'w') as lock_file:
                    fcntl.flock(lock_file, fcntl.LOCK_EX)
                    with open(error_path, "w") as f:
                        f.write(str(e))
                    # Lock is released when file is closed
            except Exception as log_error:
                # If locking fails, try direct write as fallback
                try:
                    with open(error_path, "w") as f:
                        f.write(f"{str(e)}\nAdditional error during logging: {str(log_error)}")
                except Exception:
                    # Last resort - print to stdout
                    print(f"Error processing {fname} and failed to log: {str(e)}")
    
    return completed_files


def chunk_list(lst: List, chunk_size: int) -> List[List]:
    """
    Split a list into chunks of specified size.
    
    Args:
        lst: List to split
        chunk_size: Size of each chunk
        
    Returns:
        List of chunked lists
    """
    return [lst[i:i + chunk_size] for i in range(0, len(lst), chunk_size)]


def process_files_batched(
    file_list: List[str], 
    msa_root: str, 
    save_root: str, 
    log_root: str,
    first_pdbid_to_seq: DictLike,
    seq_to_pdb_index: DictLike,
    num_workers: int,
    batch_size: Optional[int],
) -> None:
    """
    Process files in batches using parallel processing.
    
    Args:
        file_list: List of files to process
        msa_root: Root directory containing MSA files
        save_root: Root directory to save processed files
        log_root: Root directory for log files
        first_pdbid_to_seq: Dictionary mapping PDB IDs to sequences (possibly shared)
        seq_to_pdb_index: Dictionary mapping sequences to PDB indices (possibly shared)
        num_workers: Number of parallel workers to use
        batch_size: Number of files to process in each batch
    """
    if batch_size is None:
        batch_size = max(1, len(file_list) // num_workers)
    
    # Split files into batches
    batches = chunk_list(file_list, batch_size)
    print(f"Split {len(file_list)} files into {len(batches)} batches of size ~{batch_size}")
    
    # Create a partial function with fixed arguments
    batch_processor = partial(
        process_file_batch,
        msa_root=msa_root,
        save_root=save_root,
        log_root=log_root,
        first_pdbid_to_seq=first_pdbid_to_seq,
        seq_to_pdb_index=seq_to_pdb_index
    )
    
    # Track progress and timing
    start_time = time.time()
    total_processed = 0
    
    # Use ProcessPoolExecutor for parallel processing
    with concurrent.futures.ProcessPoolExecutor(max_workers=num_workers) as executor:
        # Process batches with progress tracking
        with tqdm(total=len(file_list), desc="Processing MSA files") as pbar:
            futures = []
            
            # Submit all batches to the executor
            for batch in batches:
                futures.append(executor.submit(batch_processor, batch))
            
            # Process results as they complete
            for future in concurrent.futures.as_completed(futures):
                try:
                    # Get the set of completed files
                    completed_files = future.result()
                    batch_count = len(completed_files)
                    total_processed += batch_count
                    pbar.update(batch_count)
                    
                    # Calculate and display statistics
                    elapsed = time.time() - start_time
                    files_per_second = total_processed / elapsed if elapsed > 0 else 0
                    pbar.set_postfix(
                        {"processed": total_processed, "files/sec": f"{files_per_second:.2f}"}
                    )
                except Exception as e:
                    print(f"Error processing batch: {e}")


if __name__ == "__main__":
    # Set start method to spawn to ensure compatibility with shared memory
    multiprocessing.set_start_method('spawn', force=True)
    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument("--input_msa_dir", type=str, default="./scripts/msa/data/mmcif_msa_with_taxid")
    parser.add_argument("--output_msa_dir", type=str, default="./scripts/msa/data/mmcif_msa")
    parser.add_argument("--seq_to_pdb_id", type=str, 
                        default="./scripts/msa/data/pdb_seqs/seq_to_pdb_id_entity_id.json")
    parser.add_argument("--seq_to_pdb_index", type=str, 
                        default="./scripts/msa/data/pdb_seqs/seq_to_pdb_index.json")
    parser.add_argument("--num_workers", type=int, default=multiprocessing.cpu_count(), 
                        help="Number of parallel workers to use (default: half of all CPU cores or 1)")
    parser.add_argument("--batch_size", type=int, default=None,
                        help="Number of files to process in each batch (default: adjusted automatically)")
    parser.add_argument("--shared_memory", action="store_true",
                        help="Use shared memory for dictionaries to reduce memory usage")
    args = parser.parse_args()

    msa_root = args.input_msa_dir
    save_root = args.output_msa_dir
    log_root = "./scripts/msa/data/mmcif_msa_log"
    
    # Load mapping data
    print("Loading mapping data...")
    _, first_pdbid_to_seq, seq_to_pdb_index = load_mapping_data(
        args.seq_to_pdb_id, 
        args.seq_to_pdb_index,
        use_shared_memory=args.shared_memory
    )
    print("Mapping data loaded successfully")

    os.makedirs(log_root, exist_ok=True)
    os.makedirs(save_root, exist_ok=True)

    print("Loading file names...")
    file_list = os.listdir(msa_root)
    print(f"Found {len(file_list)} MSA files to process")
    
    # Process files in batches
    process_files_batched(
        file_list=file_list,
        msa_root=msa_root,
        save_root=save_root,
        log_root=log_root,
        first_pdbid_to_seq=first_pdbid_to_seq,
        seq_to_pdb_index=seq_to_pdb_index,
        num_workers=args.num_workers,
        batch_size=args.batch_size
    )
    
    # Release all shared dictionaries if necessary
    if args.shared_memory:
        for dict_id in get_shared_dict_ids():
            release_shared_dict(dict_id)

    print("Processing complete")