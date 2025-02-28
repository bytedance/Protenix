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

import os
import mmap
import psutil
from os.path import join as opjoin
from typing import Dict, List, Tuple, Optional, Union, Any
import concurrent.futures
import multiprocessing
import time
import math

from tqdm import tqdm
from utils import (
    convert_to_shared_dict,  # To create new shared dictionaries
    release_shared_dict,     # To manually release dictionaries
    get_shared_dict_ids      # To list available dictionaries
)


def get_available_memory():
    """Get available system memory in bytes."""
    return psutil.virtual_memory().available


def estimate_memory_needs(file_size):
    """Roughly estimate memory needs for processing.
    
    Args:
        file_size (int): Size of the file in bytes
        
    Returns:
        int: Estimated memory in bytes needed for processing
    """
    # Estimate dictionary size: assume ~30% of file size for unique entries
    # plus overhead for Python dictionary (which is roughly 2-4x the raw data size)
    dict_estimate = file_size * 0.3 * 3
    
    # Some buffer for other operations
    other_memory = 100 * 1024 * 1024  # 100MB
    
    return dict_estimate + other_memory


def read_a3m(a3m_file: str) -> Tuple[List[str], List[str], int]:
    """read a3m file from output of mmseqs

    Args:
        a3m_file (str): the a3m file searched by mmseqs(colabfold search)

    Returns:
        Tuple[List[str], List[str], int]: the header, seqs of a3m files, and uniref index
    """
    heads = []
    seqs = []
    # Record the row index. The index before this index is the MSA of Uniref30 DB,
    # and the index after this index is the MSA of ColabfoldDB.
    uniref_index = 0
    with open(a3m_file, "r") as infile:
        for idx, line in enumerate(infile):
            if line.startswith(">"):
                heads.append(line)
                if idx == 0:
                    query_name = line
                elif idx > 0 and line == query_name:
                    uniref_index = idx
            else:
                seqs.append(line)
    return heads, seqs, uniref_index


def read_m8(m8_file: str) -> Dict[str, str]:
    """Read the uniref_tax.m8 file from output of mmseqs using block-wise memory mapping for efficiency.
    
    This implementation efficiently processes very large m8 files by:
    1. Using block-wise memory mapping (mmap) to only load portions of the file at a time
    2. Processing the file in multiple smaller chunks to minimize memory footprint
    3. Monitoring dictionary growth to provide relevant warnings
    
    Args:
        m8_file (str): the uniref_tax.m8 from output of mmseqs(colabfold search)

    Returns:
        Dict[str, str]: the dict mapping uniref hit_name to NCBI TaxID
    """
    # Get file size to report progress
    file_size = os.path.getsize(m8_file)
    print(f"Reading m8 file ({file_size/(1024*1024):.1f} MB)...")
    
    # Check available memory and estimate needs for the final dictionary
    available_memory = get_available_memory()
    estimated_memory = estimate_memory_needs(file_size)
    memory_ratio = estimated_memory / available_memory if available_memory > 0 else float('inf')
    
    if memory_ratio > 0.7:  # If using more than 70% of available memory
        print("⚠️ Warning: This operation may require significant memory.")
        if memory_ratio > 1.0:
            print("⚠️ Consider using the --shared_memory option for large files.")
    
    # Create a regular dictionary for fast loading
    uniref_to_ncbi_taxid = {}
    
    line_count = 0
    
    # If the file is very small, don't use mmap (overhead not worth it)
    if file_size < 10 * 1024 * 1024:  # less than 10MB
        with open(m8_file, "r") as f:
            for line in tqdm(f, desc="Reading m8 file", unit="lines"):
                line_list = line.rstrip().split("\t")
                if len(line_list) >= 3:
                    hit_name = line_list[1]
                    ncbi_taxid = line_list[2]
                    uniref_to_ncbi_taxid[hit_name] = ncbi_taxid
                line_count += 1
    else:
        # For larger files, use block-wise mmap for memory efficiency
        # Calculate block size - we'll use larger blocks (1GB) since we're only mapping one at a time
        # This provides better performance while still keeping memory usage reasonable
        block_size = 1024 * 1024 * 1024  # 1GB blocks
        num_blocks = math.ceil(file_size / block_size)
        
        # Calculate an approximate number of lines for the progress bar
        with open(m8_file, "r") as f:
            # Sample first few lines to estimate average line length
            sample_size = min(5000, file_size)
            sample = f.read(sample_size)
            line_sample_count = sample.count('\n')
            if line_sample_count > 0:
                avg_line_length = sample_size / line_sample_count
                estimated_lines = int(file_size / avg_line_length)
            else:
                estimated_lines = 1000000  # Default estimate if we can't calculate
        
        # Use a single progress bar for all blocks
        with tqdm(total=estimated_lines, desc="Reading m8 file") as pbar:
            # Process file block by block
            for block_num in range(num_blocks):
                block_offset = block_num * block_size
                block_length = min(block_size, file_size - block_offset)
                
                with open(m8_file, "r") as f:
                    # Create a memory map for just this block
                    mm = mmap.mmap(f.fileno(), length=block_length, offset=block_offset, access=mmap.ACCESS_READ)
                    
                    try:
                        # Skip partial line at the beginning if not the first block
                        if block_num > 0:
                            # Find first newline to skip partial line
                            first_newline = mm.find(b'\n')
                            if first_newline != -1:
                                mm.seek(first_newline + 1)
                        
                        # Process lines in this block
                        block_line_count = 0
                        partial_line = b''
                        
                        while True:
                            # Read a larger chunk of data for more efficient processing
                            # Increased from 1MB to 4MB chunks for better performance
                            chunk = mm.read(min(4 * 1024 * 1024, mm.size() - mm.tell()))
                            if not chunk:
                                break
                            
                            # Combine with any partial line from previous chunk
                            data = partial_line + chunk
                            
                            # Find the last newline in the chunk
                            last_newline_pos = data.rfind(b'\n')
                            
                            if last_newline_pos == -1:
                                # No newlines, store entire chunk as partial
                                partial_line = data
                                continue
                            
                            # Split the complete lines and save the partial for next iteration
                            complete_data = data[:last_newline_pos+1]
                            partial_line = data[last_newline_pos+1:]
                            
                            # Process the complete lines
                            lines = complete_data.split(b'\n')
                            for line in lines:
                                if not line:  # Skip empty lines
                                    continue
                                    
                                try:
                                    line_str = line.decode('utf-8')
                                    line_list = line_str.split('\t')
                                    
                                    # Process valid lines with at least 3 columns
                                    if len(line_list) >= 3:
                                        hit_name = line_list[1]
                                        ncbi_taxid = line_list[2]
                                        uniref_to_ncbi_taxid[hit_name] = ncbi_taxid
                                    
                                    block_line_count += 1
                                    line_count += 1
                                    
                                    # Update progress occasionally
                                    if block_line_count % 1000 == 0:
                                        pbar.update(1000)
                                except Exception as _:
                                    # Skip lines that can't be processed
                                    continue
                        
                        # Handle any remaining partial line if we're in the last block
                        if block_num == num_blocks - 1 and partial_line:
                            try:
                                line_str = partial_line.decode('utf-8')
                                line_list = line_str.split('\t')
                                
                                if len(line_list) >= 3:
                                    hit_name = line_list[1]
                                    ncbi_taxid = line_list[2]
                                    uniref_to_ncbi_taxid[hit_name] = ncbi_taxid
                                
                                line_count += 1
                            except Exception as _:
                                pass
                    
                    finally:
                        # Always close the memory map
                        mm.close()
                        
                    # Update progress bar for this block
                    if block_line_count % 1000 != 0:
                        pbar.update(block_line_count % 1000)
    
    print(f"Processed {len(uniref_to_ncbi_taxid):,} unique entries")
    
    return uniref_to_ncbi_taxid


def update_a3m(
    a3m_path: str,
    uniref_to_ncbi_taxid: Dict,
    save_root: str,
) -> str:
    """add NCBI TaxID to header if "UniRef" in header

    Args:
        a3m_path (str): the original a3m path returned by mmseqs(colabfold search)
        uniref_to_ncbi_taxid (Dict): the dict mapping uniref hit_name to NCBI TaxID
        save_root (str): the updated a3m
        
    Returns:
        str: The path of the processed a3m file
    """
    heads, seqs, uniref_index = read_a3m(a3m_path)
    fname = a3m_path.split("/")[-1]
    out_a3m_path = opjoin(save_root, fname)
    with open(out_a3m_path, "w") as ofile:
        for idx, (head, seq) in enumerate(zip(heads, seqs)):
            uniref_id = head.split("\t")[0][1:]
            ncbi_taxid = uniref_to_ncbi_taxid.get(uniref_id, None)
            if (ncbi_taxid is not None) and (idx < (uniref_index // 2)):
                if not uniref_id.startswith("UniRef100_"):
                    head = head.replace(
                        uniref_id, f"UniRef100_{uniref_id}_{ncbi_taxid}/"
                    )
                else:
                    head = head.replace(uniref_id, f"{uniref_id}_{ncbi_taxid}/")
            ofile.write(f"{head}{seq}")
    return out_a3m_path


def update_a3m_batch(batch_paths, uniref_to_ncbi_taxid, save_root):
    """Process a batch of a3m files.
    
    Args:
        batch_paths (List[str]): List of paths to a3m files to process
        uniref_to_ncbi_taxid (Dict[str, str]): Dictionary mapping UniRef IDs to NCBI TaxIDs
        save_root (str): Directory to save processed files
        
    Returns:
        List[str]: List of processed file paths
    """
    results = []
    for a3m_path in batch_paths:
        result = update_a3m(
            a3m_path=a3m_path,
            uniref_to_ncbi_taxid=uniref_to_ncbi_taxid,
            save_root=save_root
        )
        results.append(result)
    return results


def process_files(
    a3m_paths: List[str],
    uniref_to_ncbi_taxid: Union[Dict[str, str], Any],
    output_msa_dir: str,
    num_workers: Optional[int] = None,
    batch_size: Optional[int] = None
) -> None:
    """Process multiple a3m files with optimized performance.
    
    This function uses a more efficient approach for multiprocessing by using batched 
    processing to reduce the overhead of process creation and task management.
    Works with both regular dictionaries and shared dictionaries.
    
    Args:
        a3m_paths (List[str]): List of a3m file paths to process
        uniref_to_ncbi_taxid (Union[Dict[str, str], Any]): 
            Dictionary mapping UniRef IDs to NCBI TaxIDs, can be either a regular dict or a shared dict
        output_msa_dir (str): Directory to save processed files
        num_workers (int, optional): Number of worker processes. Defaults to None (uses CPU count).
        batch_size (int, optional): Size of batches for processing. Defaults to None (auto-calculated).
    """
    if num_workers is None:
        # Use a smaller number of workers to avoid excessive overhead
        num_workers = max(1, min(multiprocessing.cpu_count() - 1, 16))
    
    total_files = len(a3m_paths)
    
    if batch_size is None:
        # Calculate an optimal batch size based on number of files and workers
        # Aim for each worker to get 2-5 batches for good load balancing
        target_batches_per_worker = 3
        batch_size = max(1, math.ceil(total_files / (num_workers * target_batches_per_worker)))
    
    # Create batches
    batches = []
    for i in range(0, len(a3m_paths), batch_size):
        batch = a3m_paths[i:i + batch_size]
        batches.append(batch)
    
    # Process in single-threaded mode if we have very few files or only one worker
    if total_files < 10 or num_workers == 1:
        with tqdm(total=total_files, desc="Processing a3m files") as pbar:
            for a3m_path in a3m_paths:
                update_a3m(
                    a3m_path=a3m_path,
                    uniref_to_ncbi_taxid=uniref_to_ncbi_taxid,
                    save_root=output_msa_dir,
                )
                pbar.update(1)
        return
    
    start_time = time.time()
    
    # Use ProcessPoolExecutor for parallel processing
    with concurrent.futures.ProcessPoolExecutor(max_workers=num_workers) as executor:
        # Submit batch tasks instead of individual files
        futures = []
        for batch in batches:
            future = executor.submit(
                update_a3m_batch, 
                batch, 
                uniref_to_ncbi_taxid, 
                output_msa_dir
            )
            futures.append(future)
        
        # Track progress across all batches
        completed_files = 0
        with tqdm(total=total_files, desc="Processing a3m files") as pbar:
            for future in concurrent.futures.as_completed(futures):
                try:
                    # Each result is a list of file paths processed in the batch
                    result = future.result()
                    batch_size = len(result)
                    completed_files += batch_size
                    pbar.update(batch_size)
                except Exception as e:
                    print(f"Error processing batch: {e}")
                    # Estimate how many files might have been in this failed batch
                    avg_batch_size = total_files / len(batches)
                    pbar.update(int(avg_batch_size))
    
    end_time = time.time()
    elapsed = end_time - start_time
    print(f"Processing complete ({elapsed:.1f} seconds)")


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument("--input_msa_dir", type=str, default="./scripts/msa/data/mmcif_msa_initial")
    parser.add_argument("--output_msa_dir", type=str, default="./scripts/msa/data/mmcif_msa_with_taxid")
    parser.add_argument("--num_workers", type=int, default=None, 
                        help="Number of worker processes. Defaults to auto.")
    parser.add_argument("--batch_size", type=int, default=None,
                        help="Number of files per batch. Defaults to auto.")
    parser.add_argument("--shared_memory", action="store_true",
                        help="Use shared memory for dictionary to reduce memory usage.")
    args = parser.parse_args()
    input_msa_dir = args.input_msa_dir
    output_msa_dir = args.output_msa_dir
    os.makedirs(output_msa_dir, exist_ok=True)

    a3m_paths = os.listdir(input_msa_dir)
    a3m_paths = [opjoin(input_msa_dir, x) for x in a3m_paths if x.endswith(".a3m")]
    m8_file = f"{input_msa_dir}/uniref_tax.m8"
    
    # Read m8 file (always into a regular dictionary for performance)
    uniref_to_ncbi_taxid = read_m8(m8_file)
    
    if args.shared_memory:
        # Use shared memory dictionary
        uniref_to_ncbi_taxid = convert_to_shared_dict(uniref_to_ncbi_taxid)
    
    # Process the a3m files
    process_files(
        a3m_paths=a3m_paths,
        uniref_to_ncbi_taxid=uniref_to_ncbi_taxid,
        output_msa_dir=output_msa_dir,
        num_workers=args.num_workers,
        batch_size=args.batch_size
    )

    # Release all shared dictionaries if necessary
    if args.use_shared_memory:
        for dict_id in get_shared_dict_ids():
            release_shared_dict(dict_id)
