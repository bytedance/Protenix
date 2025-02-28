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
from os.path import join as opjoin
from typing import Dict, List, Tuple, Optional, Union, Any
import concurrent.futures
import multiprocessing
import math

from tqdm import tqdm
from utils import (
    convert_to_shared_dict,  # To create new shared dictionaries
    release_shared_dict,     # To manually release dictionaries
    get_shared_dict_ids      # To list available dictionaries
)

def process_block_binary(block_info):
    """Process a range of blocks with binary file reading for better performance
    
    Args:
        block_info (tuple): (start_block, end_block, file_path, block_size, file_size, num_blocks)
        
    Returns:
        dict: Dictionary of results from these blocks
    """
    start_block, end_block, file_path, block_size, file_size, num_blocks = block_info
    local_dict = {}
    
    # Buffer size for reading across block boundaries
    boundary_buffer_size = 8192  # 8KB should be enough for even very long lines
    
    with open(file_path, 'rb') as f:
        for block_num in range(start_block, end_block):
            # Calculate block bounds
            block_offset = block_num * block_size
            
            # Determine the start position for reading this block
            if block_num == 0:
                # First block always starts at the beginning of the file
                start_pos = 0
            else:
                # For subsequent blocks, we need to find where the first complete line starts
                # First, check if the previous block ended with a newline
                prev_block_end = block_offset - 1
                f.seek(prev_block_end)
                last_char_prev_block = f.read(1)
                
                # If the previous block ended with a newline, start at the beginning of this block
                if last_char_prev_block == b'\n':
                    start_pos = block_offset
                else:
                    # Previous block didn't end with a newline, find the first newline in this block
                    f.seek(block_offset)
                    chunk = f.read(min(boundary_buffer_size, file_size - block_offset))
                    newline_pos = chunk.find(b'\n')
                    
                    # If no newline found, this entire block is part of a line from previous block
                    if newline_pos == -1:
                        continue
                    
                    # Start after the first newline
                    start_pos = block_offset + newline_pos + 1
            
            # Calculate how much data to read from the start position
            read_length = min(block_size - (start_pos - block_offset), file_size - start_pos)
            
            # Skip if nothing to read after adjustments
            if read_length <= 0:
                continue
            
            # Read the data block
            f.seek(start_pos)
            data = f.read(read_length)
            
            # Skip if no data
            if not data:
                continue
            
            # Split into lines and process
            lines = data.split(b'\n')
            
            # Process all lines except possibly the last one if it's incomplete
            for i, line in enumerate(lines):
                # Special handling for the last line in the block (if not the last block)
                if i == len(lines) - 1 and block_num < num_blocks - 1:
                    end_pos = start_pos + len(data)
                    if end_pos < file_size and not data.endswith(b'\n'):
                        # Last line is incomplete, need to read more to complete it
                        incomplete_line = line
                        
                        # Read ahead to find the rest of the line
                        f.seek(end_pos)
                        extra_data = f.read(min(boundary_buffer_size, file_size - end_pos))
                        
                        # Find the end of the line
                        newline_pos = extra_data.find(b'\n')
                        if newline_pos != -1:
                            # Found the end of the line
                            line_remainder = extra_data[:newline_pos]
                            full_line = incomplete_line + line_remainder
                    else:
                        full_line = line
                else:
                    full_line = line

                # Process complete lines
                if full_line:  # Skip empty lines
                    try:
                        line_str = full_line.decode('utf-8')
                        line_list = line_str.split('\t')
                        hit_name = line_list[1]
                        ncbi_taxid = line_list[2]
                        local_dict[hit_name] = ncbi_taxid
                    except Exception:
                        # Skip problematic lines
                        continue
    
    return local_dict


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


def read_m8(
    m8_file: str,
    max_workers: Optional[int] = None,
    block_size_mb: int = 64
) -> Dict[str, str]:
    """Read the uniref_tax.m8 file from output of mmseqs using optimized block processing.
    
    This implementation automatically selects the best processing approach based on file size:
    1. Simple sequential processing for small to medium files
    2. Block-wise processing with multiprocessing for large files when beneficial
    
    Args:
        m8_file (str): the uniref_tax.m8 from output of mmseqs(colabfold search)
        max_workers (Optional[int]): maximum number of worker processes to use (defaults to CPU count - 1)
        block_size_mb (int): size of each processing block in MB. Do not set too small, otherwise the 
            overhead of process creation and task management will be too high and unfound bugs will be introduced.

    Returns:
        Dict[str, str]: the dict mapping uniref hit_name to NCBI TaxID
    """
    # Get file size to report progress
    file_size = os.path.getsize(m8_file)
    print(f"Reading m8 file ({file_size/(1024*1024):.1f} MB)...")
    
    # Calculate block size and number of blocks based on input parameter
    block_size = block_size_mb * 1024 * 1024  # Convert MB to bytes
    num_blocks = math.ceil(file_size / block_size)
    
    # Calculate available CPU resources
    available_cpus = multiprocessing.cpu_count() - 1 or 1  # At least 1
    
    # Determine if multiprocessing would be beneficial
    # We use multiprocessing if:
    # 1. We have more than 1 block
    # 2. We have at least 2 CPUs available
    use_multiprocessing = (
        num_blocks > 1 and 
        available_cpus > 1
    )
    
    # Set up number of workers if we're using multiprocessing
    if use_multiprocessing:
        if max_workers is None:
            # Use a reasonable number of workers based on blocks and CPUs
            max_workers = min(available_cpus, num_blocks, 16)
        else:
            # Ensure we don't use more workers than blocks or available CPUs
            max_workers = min(max_workers, available_cpus, num_blocks)
        
        # If we only have 1 worker, fall back to sequential processing
        if max_workers <= 1:
            use_multiprocessing = False
    
    uniref_to_ncbi_taxid = {}
    
    # For multiprocessing approach (large files with sufficient CPU resources)
    if use_multiprocessing:
        print(f"File is large, using multiprocessing with {max_workers} workers for {num_blocks} blocks")

        # Create batches of blocks
        batches = []
        for i in range(0, num_blocks):
            batches.append((i, i+1, m8_file, block_size, file_size, num_blocks))
        
        # Process batches with progress tracking
        with concurrent.futures.ProcessPoolExecutor(max_workers=max_workers) as executor:
            futures = [executor.submit(process_block_binary, batch) for batch in batches]
            
            with tqdm(total=len(batches), desc="Processing file blocks") as pbar:
                for future in concurrent.futures.as_completed(futures):
                    try:
                        batch_dict = future.result()
                        # Merge results
                        uniref_to_ncbi_taxid.update(batch_dict)
                        pbar.update(1)
                    except Exception as e:
                        print(f"Error processing batch: {e}")
                        pbar.update(1)

    # Sequential processing approach (for small to medium files)
    else:
        print(f"Using sequential processing for {file_size/(1024*1024):.1f} MB file")
        
        # Original single-threaded line-by-line processing
        with open(m8_file, "r") as infile:
            for line in tqdm(infile, desc="Reading m8 file", unit="lines"):
                line_list = line.rstrip().split("\t")
                hit_name = line_list[1]
                ncbi_taxid = line_list[2]
                uniref_to_ncbi_taxid[hit_name] = ncbi_taxid

    print(f"Processed {len(uniref_to_ncbi_taxid):,} unique entries")

    return uniref_to_ncbi_taxid


def update_a3m(
    a3m_path: str,
    uniref_to_ncbi_taxid: Dict[str, str],
    save_root: str,
) -> None:
    """add NCBI TaxID to header if "UniRef" in header

    Args:
        a3m_path (str): the original a3m path returned by mmseqs(colabfold search)
        uniref_to_ncbi_taxid (Dict): the dict mapping uniref hit_name to NCBI TaxID
        save_root (str): the updated a3m
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


def update_a3m_batch(batch_paths: List[str], uniref_to_ncbi_taxid: Dict[str, str], save_root: str) -> int:
    """Process a batch of a3m files.
    
    Args:
        batch_paths (List[str]): List of paths to a3m files to process
        uniref_to_ncbi_taxid (Dict[str, str]): Dictionary mapping UniRef IDs to NCBI TaxIDs
        save_root (str): Directory to save processed files
        
    Returns:
        int: Number of files processed
    """
    for a3m_path in batch_paths:
        update_a3m(
            a3m_path=a3m_path,
            uniref_to_ncbi_taxid=uniref_to_ncbi_taxid,
            save_root=save_root
        )
    return len(batch_paths)


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
    batches = [a3m_paths[i:i + batch_size] for i in range(0, len(a3m_paths), batch_size)]
    
    # Process in single-threaded mode if we have very few files or only one worker
    if total_files < 10 or num_workers == 1:
        for a3m_path in tqdm(a3m_paths, desc="Processing a3m files"):
            update_a3m(
                a3m_path=a3m_path,
                uniref_to_ncbi_taxid=uniref_to_ncbi_taxid,
                save_root=output_msa_dir,
            )
        return    
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
        with tqdm(total=total_files, desc="Processing a3m files") as pbar:
            for future in concurrent.futures.as_completed(futures):
                try:
                    # Each result is a list of file paths processed in the batch
                    batch_size = future.result()
                    pbar.update(batch_size)
                except Exception as e:
                    print(f"Error processing batch: {e}")
                    # Estimate how many files might have been in this failed batch
                    avg_batch_size = total_files / len(batches)
                    pbar.update(int(avg_batch_size))


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument("--input_msa_dir", type=str, default="./scripts/msa/data/mmcif_msa_initial")
    parser.add_argument("--output_msa_dir", type=str, default="./scripts/msa/data/mmcif_msa_with_taxid")
    parser.add_argument("--num_workers", type=int, default=None, 
                        help="Number of worker processes for a3m processing. Defaults to auto.")
    parser.add_argument("--batch_size", type=int, default=None,
                        help="Number of a3m files per batch. Defaults to auto.")
    parser.add_argument("--shared_memory", action="store_true",
                        help="Use shared memory for dictionary to reduce memory usage.")
    parser.add_argument("--mp_read_workers", type=int, default=None,
                        help="Number of worker processes for reading m8 file. Defaults to auto.")
    parser.add_argument("--block_size_mb", type=int, default=64,
                        help="Size of each processing block in MB (smaller blocks can improve parallelism).")
    args = parser.parse_args()
    
    # Set up directories
    input_msa_dir = args.input_msa_dir
    output_msa_dir = args.output_msa_dir
    os.makedirs(output_msa_dir, exist_ok=True)

    # Find input files
    a3m_paths = os.listdir(input_msa_dir)
    a3m_paths = [opjoin(input_msa_dir, x) for x in a3m_paths if x.endswith(".a3m")]
    m8_file = f"{input_msa_dir}/uniref_tax.m8"
    
    # Read m8 file with improved parameters
    print(f"Reading m8 file with block size: {args.block_size_mb}MB")
    
    uniref_to_ncbi_taxid = read_m8(
        m8_file=m8_file,
        max_workers=args.mp_read_workers,
        block_size_mb=args.block_size_mb
    )
    
    print(f"Successfully read m8 file with {len(uniref_to_ncbi_taxid):,} entries")
    
    # Convert to shared memory if needed for a3m processing
    if args.shared_memory:
        try:
            print("Converting dictionary to shared memory for a3m processing...")
            uniref_to_ncbi_taxid = convert_to_shared_dict(uniref_to_ncbi_taxid)
            print("Successfully converted dictionary to shared memory")
        except Exception as e:
            print(f"⚠️ Error converting to shared memory: {e}")
    
    # Process the a3m files
    print(f"Processing {len(a3m_paths)} a3m files...")
    process_files(
        a3m_paths=a3m_paths,
        uniref_to_ncbi_taxid=uniref_to_ncbi_taxid,
        output_msa_dir=output_msa_dir,
        num_workers=args.num_workers,
        batch_size=args.batch_size
    )
    
    # Release all shared dictionaries if necessary
    if args.shared_memory:
        for dict_id in get_shared_dict_ids():
            try:
                release_shared_dict(dict_id)
            except Exception as e:
                print(f"Warning: Failed to release shared dict {dict_id}: {e}")

    print("Processing complete")
