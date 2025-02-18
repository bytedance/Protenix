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

import argparse
import functools
import logging
import multiprocessing
import pickle
import subprocess as sp
from pathlib import Path
from typing import Optional, Union

import gemmi
import numpy as np
import rdkit
import tqdm
from biotite.structure.io import pdbx
from pdbeccdutils.core import ccd_reader


def download_ccd_cif(output_path: Path):
    """
    Download the CCD CIF file from rcsb.org.

    Args:
        output_path (Path): The output path for saving the downloaded CCD CIF file.
    """
    output_path.parent.mkdir(parents=True, exist_ok=True)

    logging.info("Downloading CCD CIF file from rcsb.org ...")

    output_cif_gz = output_path / "components.cif.gz"
    if output_cif_gz.exists():
        logging.info("Remove old zipped CCD CIF file: %s", output_cif_gz)
        output_cif_gz.unlink()

    output_cif = output_cif_gz.with_suffix("")
    if output_cif.exists():
        logging.info("Remove old CCD CIF file: %s", output_cif)
        output_cif.unlink()

    sp.run(
        f"wget https://files.wwpdb.org/pub/pdb/data/monomers/components.cif.gz -P {output_path}",
        shell=True,
        check=True,
    )

    sp.run(
        f"gunzip {output_cif_gz}",
        shell=True,
        check=True,
    )

    # remove .gz file
    output_cif_gz.unlink(missing_ok=True)

    logging.info("Download CCD CIF file successfully.")


@functools.lru_cache
def gemmi_load_ccd_cif(ccd_cif: Union[Path, str]) -> gemmi.cif.Document:
    """
    Load CCD components file by gemmi

    ccd_cif (Union[Path, str]): The path to the CCD CIF file.

    Returns:
        Document: gemmi ccd components file
    """
    return gemmi.cif.read(str(ccd_cif))


@functools.lru_cache
def biotite_load_ccd_cif(ccd_cif: Union[Path, str]) -> pdbx.CIFFile:
    """
    Load CCD components file by biotite

    Args:
        ccd_cif (Union[Path, str]): The path to the CCD CIF file.

    Returns:
        pdbx.CIFFile: ccd components file
    """
    return pdbx.CIFFile.read(str(ccd_cif))


def _get_component_rdkit_mol_processing(
    ccd_code_and_cif_file: tuple[str, Path]
) -> Optional[rdkit.Chem.Mol]:
    """
    Get rdkit mol by PDBeCCDUtils
    https://github.com/PDBeurope/ccdutils

    Args:
        ccd_code (str): ccd code
        ccd_cif_file (Path): The path to the CCD CIF file.

    Returns
        rdkit.Chem.Mol: rdkit mol with ref coord
    """
    ccd_code, ccd_cif_file = ccd_code_and_cif_file
    ccd_cif = gemmi_load_ccd_cif(ccd_cif_file)
    try:
        ccd_block = ccd_cif[ccd_code]
    except KeyError:
        return None
    ccd_reader_result = ccd_reader._parse_pdb_mmcif(ccd_block, sanitize=True)
    mol = ccd_reader_result.component.mol

    # atom name from ccd, reading by pdbeccdutils
    # copy atom name for pickle https://github.com/rdkit/rdkit/issues/2470
    mol.atom_map = {atom.GetProp("name"): atom.GetIdx() for atom in mol.GetAtoms()}

    mol.name = ccd_code
    mol.sanitized = ccd_reader_result.sanitized
    mol.ref_conf_id = 0  # first conf is ideal conf.
    mol.ref_conf_type = "idea"

    num_atom = mol.GetNumAtoms()
    if num_atom == 0:  # eg: UNL without atom
        return mol

    # make ref_mask, ref_mask is True if ideal coord is valid
    atoms = ccd_block.find(
        "_chem_comp_atom.", ["atom_id", "model_Cartn_x", "pdbx_model_Cartn_x_ideal"]
    )
    assert num_atom == len(atoms)
    ref_mask = np.zeros(num_atom, dtype=bool)
    for row in atoms:
        atom_id = gemmi.cif.as_string(row["_chem_comp_atom.atom_id"])
        atom_idx = mol.atom_map[atom_id]
        x_ideal = row["_chem_comp_atom.pdbx_model_Cartn_x_ideal"]
        ref_mask[atom_idx] = x_ideal != "?"
    mol.ref_mask = ref_mask

    if mol.sanitized == False:
        return mol
    options = rdkit.Chem.AllChem.ETKDGv3()
    options.clearConfs = False
    try:
        conf_id = rdkit.Chem.AllChem.EmbedMolecule(mol, options)
        mol.ref_conf_id = conf_id
        mol.ref_conf_type = "rdkit"
        mol.ref_mask[:] = True
    except Exception:
        logging.warning(
            "Warning: fail to generate conf for %s, use idea conf", ccd_code
        )  # sanitization issue here
    return mol


def precompute_ccd_mol(ccd_cif: Path, output_pkl: Path, num_cpu: int = 1):
    """
    Precompute the CCD CIF file.

    Args:
        cif_file (Path): The path to the CCD CIF file.
        output_pkl (Path): The output path for saving the precomputed CCD CIF file.
        num_cpu (int): The number of CPUs to use for parallel processing.
    """
    # preprocessing all ccd components in _components_file at first time run.
    gemmi_load_ccd_cif(ccd_cif)

    mols = {}

    biotite_ccd_cif = biotite_load_ccd_cif(ccd_cif)
    ccd_codes = list(biotite_ccd_cif.keys())

    tasks = list(zip(ccd_codes, [ccd_cif] * len(ccd_codes)))

    with multiprocessing.Pool(num_cpu) as pool:
        for mol in tqdm.tqdm(
            pool.imap_unordered(
                _get_component_rdkit_mol_processing,
                tasks,
            ),
            smoothing=0,
            total=len(ccd_codes),
        ):
            if mol is None:
                continue
            mols[mol.name] = mol

    # success rate
    n_ccd = len(ccd_codes)
    logging.info(
        "success rate: %.2f%% (%d/%d)", len(mols) / n_ccd * 100, len(mols), n_ccd
    )

    # sanitized rate
    sanitized_num = sum([mol.sanitized for mol in mols.values()])
    logging.info(
        "sanitized rate: %.2f%% (%d/%d)",
        sanitized_num / n_ccd * 100,
        sanitized_num,
        n_ccd,
    )

    # rdkit conf rate
    rdkit_conf_num = sum([mol.ref_conf_type == "rdkit" for mol in mols.values()])
    logging.info(
        "rdkit conf rate: %.2f%% (%d/%d)",
        rdkit_conf_num / n_ccd * 100,
        rdkit_conf_num,
        n_ccd,
    )

    with open(output_pkl, "wb") as f:
        pickle.dump(mols, f)
    logging.info("save rdkit mol to %s", output_pkl)

    ccd_list_txt = ccd_cif.with_suffix(".txt")
    with open(ccd_list_txt, "w") as f:
        f.write("\n".join(mols.keys()))


def run_update_ccd_cache(
    ccd_cache_dir: Path, num_cpu: int = 1, disable_download: bool = False
):
    """
    Updates the CCD (Chemical Component Dictionary) cache by downloading the latest
    CCD CIF file and precomputing RDKit molecule objects.

    Args:
        ccd_cache_dir (Path): The directory where the CCD cache files are stored.
        num_cpu (int, optional): The number of CPU cores to use for precomputing RDKit molecules.
                                 Defaults to 1.
        disable_download (bool, optional): If True, skips downloading the CCD CIF file.
                                           Defaults to False.
    """

    if not disable_download:
        download_ccd_cif(output_path=ccd_cache_dir)

    ccd_cif = ccd_cache_dir / "components.cif"
    ccd_rdkit_mol_pkl = ccd_cache_dir / "components.cif.rdkit_mol.pkl"
    precompute_ccd_mol(ccd_cif, ccd_rdkit_mol_pkl, num_cpu=num_cpu)


if __name__ == "__main__":

    current_file_path = Path(__file__)
    current_directory = current_file_path.parent
    code_directory = current_directory.parent
    releases_data_ccd_directory = code_directory / "release_data" / "ccd_cache"

    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-c",
        "--ccd_cache_dir",
        type=Path,
        default=releases_data_ccd_directory,
        help='Path to the CCD cache directory. Defaults to "release_data/ccd_cache" under the code directory.',
    )
    parser.add_argument(
        "-n",
        "--n_cpu",
        type=int,
        default=1,
        help="Number of worker processes to use. Defaults to 1.",
    )

    parser.add_argument(
        "-d",
        "--disable_download",
        action="store_true",
        help="Whether to disable downloading the CCD CIF file. Defaults to False.",
    )

    args = parser.parse_args()

    run_update_ccd_cache(
        ccd_cache_dir=args.ccd_cache_dir,
        num_cpu=args.n_cpu,
        disable_download=args.disable_download,
    )
