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

import gzip
import json
import pickle
import time
from typing import Any, TYPE_CHECKING, Union

# Only import lmdb for type checking (no runtime cost)
if TYPE_CHECKING:
    import lmdb

import hashlib
import os
from pathlib import Path

import pandas as pd
import torch
from biotite import structure
from biotite.structure.io import pdbx

from protenix.data.constants import mmcif_restype_1to3
from protenix.utils.logger import get_logger
from protenix.utils.torch_utils import map_values_to_list, to_device

logger = get_logger(__name__)

PANDAS_NA_VALUES = [
    "",
    "#N/A",
    "#N/A N/A",
    "#NA",
    "-1.#IND",
    "-1.#QNAN",
    "-NaN",
    "-nan",
    "1.#IND",
    "1.#QNAN",
    "<NA>",
    "N/A",
    # "NA",
    "NULL",
    "NaN",
    "n/a",
    "nan",
    "null",
]

_INDICES_CSV_CACHE: dict[str, pd.DataFrame] = {}
_JSON_FILE_CACHE: dict[str, Any] = {}


def read_indices_csv(
    csv: Union[str, Path], load_indice_mode: str = "full"
) -> pd.DataFrame:
    """
    Read a csv file without the content changing.

    Args:
        csv (Union[str, Path]): A csv file path.
        load_indice_mode (str): Mode to load indices. Defaults to "full".

    Returns:
        pd.DataFrame : A pandas DataFrame.
    """
    csv_str = str(csv)

    if csv_str in _INDICES_CSV_CACHE:
        df_cached = _INDICES_CSV_CACHE[csv_str]
        logger.info(
            "[Protenix IO] read_indices_csv cache hit: path=%s, shape=%s",
            csv_str,
            tuple(df_cached.shape),
        )
        # Always return a copy so per-dataset filtering does not mutate the cache.
        return df_cached.copy(deep=True)

    t0 = time.time()
    if csv_str.endswith(".parquet"):
        if load_indice_mode == "simple":
            df = pd.read_parquet(
                csv_str,
                columns=["pdb_id", "type", "chain_1_id", "chain_2_id", "cluster_id"],
            )
        else:
            df = pd.read_parquet(csv_str)
        df = df.astype(str)
    else:
        df = pd.read_csv(
            csv_str, na_values=PANDAS_NA_VALUES, keep_default_na=False, dtype=str
        )

    def _normalize_entity_id(entity_id):
        if pd.isna(entity_id):
            return None
        try:
            return str(int(float(entity_id)))
        except (ValueError, TypeError):
            return None

    if "entity_1_id" in df.columns:
        df["entity_1_id"] = df["entity_1_id"].apply(_normalize_entity_id)
    if "entity_2_id" in df.columns:
        df["entity_2_id"] = df["entity_2_id"].apply(_normalize_entity_id)

    _INDICES_CSV_CACHE[csv_str] = df
    logger.info(
        "[Protenix IO] read_indices_csv finished in %.3fs: path=%s, shape=%s",
        time.time() - t0,
        csv_str,
        tuple(df.shape),
    )
    # Return a copy so callers can safely modify.
    return df.copy(deep=True)


class LMDBDict:
    """
    A lightweight dict-like wrapper around an LMDB database using SHA1-hashed keys.
    Values are read lazily from disk on access.

    Only record path and configuration at init time; do not open the
    Environment yet. Each process will open its own read-only env on the
    first access (safe under fork-based multiprocessing).

    Args:
        lmdb_path: Path to the LMDB database.
        lock (bool, optional): Whether to use locking. Defaults to False.
    """

    def __init__(self, lmdb_path, lock=False):
        self.path = str(lmdb_path)
        self.lock = lock
        # Per-process cache: pid -> lmdb.Environment
        self._env_by_pid: dict[int, Any] = {}  # Avoid runtime dependency on lmdb.Env

    def _get_env(self) -> "lmdb.Environment":
        """
        Maintain one read-only env per process to avoid sharing an Environment
        across forked processes.
        """
        pid = os.getpid()
        env = self._env_by_pid.get(pid)
        if env is None:
            # Lazy import: only import lmdb when actually needed
            try:
                import lmdb
            except ImportError as e:
                raise ImportError(
                    "LMDB is required to use LMDBDict. Please install it via: pip install lmdb"
                ) from e

            env = lmdb.open(
                self.path,
                subdir=False,
                readonly=True,
                lock=self.lock,
                readahead=False,
                meminit=False,
            )
            self._env_by_pid[pid] = env
        return env

    def _hash_key(self, key):
        """Compute SHA1 hash for a key, consistent with the writer script."""
        return hashlib.sha1(str(key).encode("utf-8")).hexdigest().encode("utf-8")

    def __getitem__(self, key):
        hashed_key = self._hash_key(key)
        env = self._get_env()
        with env.begin() as txn:
            value_bytes = txn.get(hashed_key)
            if value_bytes is None:
                raise KeyError(f"Key {key} not found in LMDB")
            # Assume the stored value is a UTF-8 string; fall back to raw bytes
            # to support arbitrary binary payloads.
            try:
                return value_bytes.decode("utf-8")
            except UnicodeDecodeError:
                return value_bytes

    def get(self, key, default=None):
        try:
            return self.__getitem__(key)
        except KeyError:
            return default

    def __contains__(self, key):
        hashed_key = self._hash_key(key)
        env = self._get_env()
        with env.begin() as txn:
            return txn.get(hashed_key) is not None

    def __len__(self):
        env = self._get_env()
        with env.begin() as txn:
            return txn.stat()["entries"]

    def keys(self):
        # Warning: do not convert keys() of a huge LMDB into a list().
        env = self._get_env()
        with env.begin() as txn:
            cursor = txn.cursor()
            for k, _ in cursor:
                # Note: yields hashed keys, original keys cannot be recovered here.
                yield k

    def close(self):
        for env in self._env_by_pid.values():
            try:
                env.close()
            except Exception:
                pass
        self._env_by_pid.clear()


def load_json_cached(path: Union[str, Path]) -> Any:
    """
    Load a JSON/Parquet/LMDB file with a simple in-process cache.

    - For JSON: Loads entire file into RAM (returns dict).
    - For Parquet: Loads into DataFrame (returns DataFrame).
    - For LMDB: Opens connection handle (returns LMDBDict wrapper).
    """
    path_str = str(Path(path))

    # 1. Cache JSON results per process; LMDB always returns a fresh wrapper.
    if path_str.endswith(".json") and path_str in _JSON_FILE_CACHE:
        logger.info("[Protenix IO] load_json_cached cache hit: path=%s", path_str)
        return _JSON_FILE_CACHE[path_str]

    t0 = time.time()

    # 2. Load according to suffix
    if path_str.endswith(".lmdb"):
        # LMDB branch: construct a new LMDBDict instance (not cached in _JSON_FILE_CACHE).
        data = LMDBDict(path_str)

    elif path_str.endswith(".json"):
        # JSON branch: load full JSON into memory and cache it.
        with open(path_str, "r") as f:
            try:
                import orjson

                data = orjson.loads(f.read())
            except ImportError:
                data = json.load(f)
        _JSON_FILE_CACHE[path_str] = data
    elif path_str == ".":
        # Empty JSON file
        data = {}
    else:
        raise ValueError(f"Unsupported file format: {path_str}")

    logger.info(
        "[Protenix IO] load_json_cached finished in %.3fs: path=%s",
        time.time() - t0,
        path_str,
    )
    return data


def load_gzip_pickle(pkl: Union[str, Path]) -> Any:
    """
    Load a gzip pickle file.

    Args:
        pkl (Union[str, Path]): A gzip pickle file path.

    Returns:
        Any: The loaded data.
    """
    with gzip.open(pkl, "rb") as f:
        data = pickle.load(f)
    return data


def dump_gzip_pickle(data: Any, pkl: Union[str, Path]):
    """
    Dump a gzip pickle file.

    Args:
        data (Any): The data to be dumped.
        pkl (Union[str, Path]): A gzip pickle file path.
    """
    with gzip.open(pkl, "wb") as f:
        pickle.dump(data, f)


def pdb_to_cif(pdb_file: str, cif_file: str, sequence=None):
    """
    Convert monomer protein pdb to CIF, and fill "entity_poly" and "entity_poly_seq"category.

    pdb_file (str): PDB file path.
    cif_file (str): output CIF file path.
    """
    entity_poly = pdbx.CIFCategory({"entity_id": "0", "type": "polypeptide(L)"})
    block_dict = {"entity_poly": entity_poly}

    if sequence is not None:
        entity_id = ["0"] * len(sequence)
        num = []
        mon_id = []
        for idx, i in enumerate(sequence):
            num.append(str(idx + 1))
            mon_id.append(mmcif_restype_1to3.get(i, "UNK"))
            assert (
                i in mmcif_restype_1to3 or i == "X"
            ), f"{i} is not in mmcif_restype_1to3"

        block_dict["entity_poly_seq"] = pdbx.CIFCategory(
            {
                "entity_id": entity_id,
                "num": num,
                "mon_id": mon_id,
            }
        )

    block = pdbx.CIFBlock(block_dict)
    pdb_file = Path(pdb_file)
    block_name = pdb_file.stem
    cif = pdbx.CIFFile({block_name: block})
    atom_array = structure.io.load_structure(pdb_file, extra_fields=["b_factor"])
    pdbx.set_structure(cif, atom_array)
    cif.write(cif_file)


class FloatEncoder(json.JSONEncoder):
    def __init__(self, *args, **kwargs):
        self.precision = kwargs.pop("precision", 2)
        super(FloatEncoder, self).__init__(*args, **kwargs)

    def encode(self, obj):
        def float_converter(o):
            if isinstance(o, float):
                return format(o, f".{self.precision}f")
            if isinstance(o, list):
                return [float_converter(i) for i in o]
            if isinstance(o, dict):
                return {k: float_converter(v) for k, v in o.items()}
            return o

        return super(FloatEncoder, self).encode(float_converter(obj))


def save_json(data, output_fpath, indent=4):
    data_json = data.copy()
    data_json = map_values_to_list(data_json)
    with open(output_fpath, "w") as f:
        if indent is not None:
            json.dump(data_json, f, indent=indent)
        else:
            json.dump(data_json, f)


def save_tensor(data, output_fpath):
    torch.save(to_device(data, device=None), output_fpath)
