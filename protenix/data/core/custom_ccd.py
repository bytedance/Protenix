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

import tempfile
from pathlib import Path
from typing import Any, Mapping

import biotite
import gemmi
import numpy as np
import rdkit
from biotite.structure import AtomArray
from biotite.structure.io import pdbx
from pdbeccdutils.core import ccd_reader
from rdkit import Chem
from rdkit.Chem import AllChem

from protenix.data.core import ccd


class CCDProvider:
    """Per-job CCD overlay with fallback to the global Protenix CCD cache."""

    def __init__(self, user_ccd_path: Path | None = None) -> None:
        self.user_ccd_path = user_ccd_path
        self._user_biotite_cif: pdbx.CIFFile | None = None
        self._user_gemmi_cif: gemmi.cif.Document | None = None
        self._custom_mols: dict[str, Chem.Mol] = {}
        self._tmpdir: tempfile.TemporaryDirectory | None = None
        if user_ccd_path is not None:
            self._load_user_ccd(user_ccd_path)

    @classmethod
    def from_job(
        cls,
        job: Mapping[str, Any],
        base_dir: str | Path | None = None,
    ) -> "CCDProvider":
        """Build a provider from optional job-level userCCD/userCCDPath fields."""
        has_user_ccd = "userCCD" in job and job["userCCD"] is not None
        has_user_ccd_path = "userCCDPath" in job and job["userCCDPath"] is not None
        if has_user_ccd and has_user_ccd_path:
            raise ValueError('Only one of "userCCD" and "userCCDPath" may be set.')
        if has_user_ccd:
            user_ccd = str(job["userCCD"])
            if not user_ccd.strip():
                raise ValueError('"userCCD" can not be empty.')
            tmpdir = tempfile.TemporaryDirectory(prefix="protenix_user_ccd_")
            path = Path(tmpdir.name) / "user_components.cif"
            path.write_text(user_ccd, encoding="utf-8")
            provider = cls(path)
            provider._tmpdir = tmpdir
            return provider
        if has_user_ccd_path:
            user_ccd_path = str(job["userCCDPath"])
            if not user_ccd_path.strip():
                raise ValueError('"userCCDPath" can not be empty.')
            path = Path(user_ccd_path).expanduser()
            if not path.is_absolute() and base_dir is not None:
                path = Path(base_dir) / path
            if not path.exists() or not path.is_file():
                raise FileNotFoundError(f'userCCDPath does not exist: "{path}"')
            return cls(path)
        return cls()

    def has_user_components(self) -> bool:
        return self._user_biotite_cif is not None

    @property
    def user_codes(self) -> set[str]:
        if self._user_biotite_cif is None:
            return set()
        return set(self._user_biotite_cif.keys())

    def get_custom_mols(self) -> dict[str, Chem.Mol]:
        return self._custom_mols

    def get_component_atom_array(
        self,
        ccd_code: str,
        keep_leaving_atoms: bool = False,
        keep_hydrogens: bool = False,
    ) -> AtomArray | None:
        if ccd_code not in self.user_codes:
            return ccd.get_component_atom_array(
                ccd_code,
                keep_leaving_atoms=keep_leaving_atoms,
                keep_hydrogens=keep_hydrogens,
            )

        assert self._user_biotite_cif is not None
        try:
            comp = pdbx.get_component(
                self._user_biotite_cif,
                data_block=ccd_code,
                use_ideal_coord=True,
                allow_missing_coord=True,
            )
        except biotite.InvalidFileError as exc:
            raise ValueError(
                f"Can not parse user CCD component {ccd_code}: {exc}"
            ) from exc

        atom_category = self._user_biotite_cif[ccd_code]["chem_comp_atom"]
        try:
            leaving_atom_flag = atom_category["pdbx_leaving_atom_flag"].as_array()
        except KeyError:
            leaving_atom_flag = np.array(["N"] * len(comp))
        comp.set_annotation("leaving_atom_flag", leaving_atom_flag == "Y")

        for atom_id in ["alt_atom_id", "pdbx_component_atom_id"]:
            try:
                comp.set_annotation(atom_id, atom_category[atom_id].as_array())
            except KeyError:
                comp.set_annotation(atom_id, comp.atom_name.copy())

        if not keep_leaving_atoms:
            comp = comp[~comp.leaving_atom_flag]
        if not keep_hydrogens:
            comp = comp[~np.isin(comp.element, ["H", "D"])]

        comp.central_to_leaving_groups = ccd._map_central_to_leaving_groups(comp)
        return comp

    def get_mol_type(self, ccd_code: str) -> str:
        if ccd_code not in self.user_codes:
            return ccd.get_mol_type(ccd_code)

        assert self._user_biotite_cif is not None
        link_type = (
            self._user_biotite_cif[ccd_code]["chem_comp"]["type"].as_item().upper()
        )
        if "PEPTIDE" in link_type and link_type != "PEPTIDE-LIKE":
            return "protein"
        if "DNA" in link_type:
            return "dna"
        if "RNA" in link_type:
            return "rna"
        return "ligand"

    def get_one_letter_code(self, ccd_code: str) -> str | None:
        if ccd_code not in self.user_codes:
            return ccd.get_one_letter_code(ccd_code)

        assert self._user_biotite_cif is not None
        one = self._user_biotite_cif[ccd_code]["chem_comp"]["one_letter_code"].as_item()
        if one == "?":
            return None
        return one

    def get_component_rdkit_mol(self, ccd_code: str) -> Chem.Mol | None:
        if ccd_code not in self.user_codes:
            return ccd.get_component_rdkit_mol(ccd_code)
        if ccd_code not in self._custom_mols:
            self._custom_mols[ccd_code] = self._build_user_rdkit_mol(ccd_code)
        return self._custom_mols[ccd_code]

    def get_ccd_ref_info(
        self,
        ccd_code: str,
        return_perm: bool = True,
        return_atomic_number: bool = False,
    ) -> dict[str, Any]:
        if ccd_code not in self.user_codes:
            return ccd.get_ccd_ref_info(
                ccd_code,
                return_perm=return_perm,
                return_atomic_number=return_atomic_number,
            )
        mol = self.get_component_rdkit_mol(ccd_code)
        if mol is None:
            return {}
        return ccd.get_ccd_ref_info(
            ccd_code,
            return_perm=return_perm,
            ccd_mols=((ccd_code, mol),),
            return_atomic_number=return_atomic_number,
        )

    def _load_user_ccd(self, user_ccd_path: Path) -> None:
        try:
            self._user_biotite_cif = pdbx.CIFFile.read(str(user_ccd_path))
            self._user_gemmi_cif = gemmi.cif.read(str(user_ccd_path))
        except Exception as exc:
            raise ValueError(
                f"Failed to parse user CCD file {user_ccd_path}: {exc}"
            ) from exc
        if not list(self._user_biotite_cif.keys()):
            raise ValueError(
                f"User CCD file {user_ccd_path} does not contain components."
            )
        for code in self._user_biotite_cif.keys():
            self._validate_user_component(code)

    def _validate_user_component(self, ccd_code: str) -> None:
        assert self._user_biotite_cif is not None
        block = self._user_biotite_cif[ccd_code]
        for category in ["chem_comp", "chem_comp_atom", "chem_comp_bond"]:
            try:
                block[category]
            except KeyError as exc:
                raise ValueError(
                    f'User CCD component "{ccd_code}" is missing required '
                    f'category "{category}".'
                ) from exc
        atom_category = block["chem_comp_atom"]
        atom_names = atom_category["atom_id"].as_array()
        if len(set(atom_names)) != len(atom_names):
            raise ValueError(
                f'User CCD component "{ccd_code}" has duplicate atom names.'
            )

    def _build_user_rdkit_mol(self, ccd_code: str) -> Chem.Mol:
        assert self._user_gemmi_cif is not None
        try:
            ccd_block = self._user_gemmi_cif[ccd_code]
        except KeyError as exc:
            raise ValueError(f'User CCD component "{ccd_code}" not found.') from exc

        try:
            ccd_reader_result = ccd_reader._parse_pdb_mmcif(ccd_block, sanitize=True)
        except Exception as exc:
            raise ValueError(
                "Failed to build RDKit molecule for user CCD component "
                f'"{ccd_code}": {exc}'
            ) from exc
        mol = ccd_reader_result.component.mol
        mol.atom_map = {atom.GetProp("name"): atom.GetIdx() for atom in mol.GetAtoms()}
        mol.name = ccd_code
        mol.sanitized = ccd_reader_result.sanitized
        mol.ref_conf_id = 0
        mol.ref_conf_type = "ideal"

        num_atom = mol.GetNumAtoms()
        if num_atom == 0:
            mol.ref_mask = np.zeros(0, dtype=bool)
            return mol

        mol.ref_mask = self._build_ref_mask(ccd_block, mol)

        if not mol.sanitized:
            return mol
        options = AllChem.ETKDGv3()
        options.clearConfs = False
        try:
            conf_id = AllChem.EmbedMolecule(mol, options)
            if conf_id >= 0:
                mol.ref_conf_id = conf_id
                mol.ref_conf_type = "rdkit"
                mol.ref_mask[:] = True
        except Exception:
            pass
        return mol

    @staticmethod
    def _build_ref_mask(ccd_block: gemmi.cif.Block, mol: rdkit.Chem.Mol) -> np.ndarray:
        ref_mask = np.ones(mol.GetNumAtoms(), dtype=bool)
        atoms = ccd_block.find(
            "_chem_comp_atom.", ["atom_id", "pdbx_model_Cartn_x_ideal"]
        )
        if len(atoms) != mol.GetNumAtoms():
            return ref_mask
        ref_mask[:] = False
        for row in atoms:
            atom_id = gemmi.cif.as_string(row["_chem_comp_atom.atom_id"])
            atom_idx = mol.atom_map[atom_id]
            x_ideal = row["_chem_comp_atom.pdbx_model_Cartn_x_ideal"]
            ref_mask[atom_idx] = x_ideal != "?"
        return ref_mask


DEFAULT_CCD_PROVIDER = CCDProvider()
