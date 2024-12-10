# Copyright 2024 ByteDance and/or its affiliates.
#
# Licensed under Creative Commons Attribution-NonCommercial 4.0
# International License (the "License");  you may not use this file  except
# in compliance with the License. You may obtain a copy of the License at
#
#     http://creativecommons.org/licenses/by-nc/4.0/
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

import biotite.structure as struc
import numpy as np
from biotite.structure import AtomArray

from protenix.data.constants import CRYSTALLIZATION_AIDS


class Filter(object):
    """
    Ref: AlphaFold3 SI Chapter 2.5.4
    """

    @staticmethod
    def remove_hydrogens(atom_array: AtomArray) -> AtomArray:
        """remove hydrogens and deuteriums"""
        return atom_array[~np.isin(atom_array.element, ["H", "D"])]

    @staticmethod
    def remove_water(atom_array: AtomArray) -> AtomArray:
        """remove water (HOH) and deuterated water (DOD)"""
        return atom_array[~np.isin(atom_array.res_name, ["HOH", "DOD"])]

    @staticmethod
    def remove_element_X(atom_array: AtomArray) -> AtomArray:
        """
        remove element X
        following residues have element X:
        - UNX: unknown one atom or ion
        - UNL: unknown ligand, some atoms are marked as X
        - ASX: ASP/ASN ambiguous, two ambiguous atoms are marked as X, 6 entries in the PDB
        - GLX: GLU/GLN ambiguous, two ambiguous atoms are marked as X, 5 entries in the PDB
        """
        X_mask = np.zeros(len(atom_array), dtype=bool)
        starts = struc.get_residue_starts(atom_array, add_exclusive_stop=True)
        for start, stop in zip(starts[:-1], starts[1:]):
            res_name = atom_array.res_name[start]
            if res_name in ["UNX", "UNL"]:
                X_mask[start:stop] = True
        atom_array = atom_array[~X_mask]

        # map ASX to ASP, as ASP is more symmetric than ASN
        mask = atom_array.res_name == "ASX"
        atom_array.res_name[mask] = "ASP"
        atom_array.atom_name[mask & (atom_array.atom_name == "XD1")] = "OD1"
        atom_array.atom_name[mask & (atom_array.atom_name == "XD2")] = "OD2"
        atom_array.element[mask & (atom_array.element == "X")] = "O"

        # map GLX to GLU, as GLU is more symmetric than GLN
        mask = atom_array.res_name == "GLX"
        atom_array.res_name[mask] = "GLU"
        atom_array.atom_name[mask & (atom_array.atom_name == "XE1")] = "OE1"
        atom_array.atom_name[mask & (atom_array.atom_name == "XE2")] = "OE2"
        atom_array.element[mask & (atom_array.element == "X")] = "O"
        return atom_array

    @staticmethod
    def remove_crystallization_aids(
        atom_array: AtomArray, entity_poly_type: dict
    ) -> AtomArray:
        """remove crystallization aids, eg: SO4, GOL, etc.

        Only remove crystallization aids if the chain is not polymer.

        Ref: AlphaFold3 SI Chapter 2.5.4
        """
        non_aids_mask = ~np.isin(atom_array.res_name, CRYSTALLIZATION_AIDS)
        poly_mask = np.isin(atom_array.label_entity_id, list(entity_poly_type.keys()))
        return atom_array[poly_mask | non_aids_mask]
