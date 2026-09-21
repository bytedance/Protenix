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
import unittest
from pathlib import Path

import numpy as np
from biotite.structure import AtomArray

from protenix.data.inference.structure_scoring import (
    build_chain_id_map,
    build_source_coord_map,
    collect_structure_files,
    copy_label_asym_ids_to_entity_ids,
    map_source_coords_to_internal,
)


def _atom_array(chain_ids, res_ids, atom_names, res_names, coords):
    atom_array = AtomArray(len(chain_ids))
    atom_array.chain_id = np.array(chain_ids)
    atom_array.res_id = np.array(res_ids)
    atom_array.atom_name = np.array(atom_names)
    atom_array.res_name = np.array(res_names)
    atom_array.coord = np.asarray(coords, dtype=np.float32)
    return atom_array


class TestStructureScoring(unittest.TestCase):
    def test_copy_label_asym_ids_to_entity_ids(self):
        input_json = {
            "sequences": [
                {
                    "proteinChain": {
                        "sequence": "AA",
                        "count": 2,
                        "label_asym_id": ["A", "B"],
                    }
                },
                {"ligand": {"ligand": "CCD_ATP", "count": 1, "label_asym_id": ["C"]}},
            ],
            "name": "sample",
        }

        self.assertTrue(copy_label_asym_ids_to_entity_ids(input_json))
        self.assertEqual(input_json["sequences"][0]["proteinChain"]["id"], ["A", "B"])
        self.assertEqual(input_json["sequences"][1]["ligand"]["id"], ["C"])

    def test_copy_label_asym_ids_rejects_duplicates(self):
        input_json = {
            "sequences": [
                {
                    "proteinChain": {
                        "sequence": "AA",
                        "count": 2,
                        "label_asym_id": ["A", "A"],
                    }
                },
            ],
            "name": "sample",
        }

        self.assertFalse(copy_label_asym_ids_to_entity_ids(input_json))
        self.assertNotIn("id", input_json["sequences"][0]["proteinChain"])

    def test_collect_structure_files_filters_supported_suffixes(self):
        with tempfile.TemporaryDirectory() as tmp_dir:
            root = Path(tmp_dir)
            (root / "a.pdb").write_text("", encoding="utf-8")
            (root / "b.cif").write_text("", encoding="utf-8")
            (root / "ignore.txt").write_text("", encoding="utf-8")
            nested = root / "nested"
            nested.mkdir()
            (nested / "c.pdb").write_text("", encoding="utf-8")

            files = collect_structure_files(root)
            self.assertEqual([path.name for path in files], ["a.pdb", "b.cif"])

            recursive_files = collect_structure_files(root, recursive=True)
            self.assertEqual(
                [path.name for path in recursive_files],
                ["a.pdb", "b.cif", "c.pdb"],
            )

    def test_build_chain_id_map_prefers_identity_over_order(self):
        source = _atom_array(
            chain_ids=["B", "A"],
            res_ids=[1, 1],
            atom_names=["N", "N"],
            res_names=["GLY", "ALA"],
            coords=[[2, 0, 0], [1, 0, 0]],
        )
        internal = _atom_array(
            chain_ids=["A", "B"],
            res_ids=[1, 1],
            atom_names=["N", "N"],
            res_names=["ALA", "GLY"],
            coords=[[1, 0, 0], [2, 0, 0]],
        )

        chain_id_map = build_chain_id_map(source, internal)

        self.assertEqual(chain_id_map["internal_to_source"], {"A": "A", "B": "B"})

    def test_build_chain_id_map_can_fall_back_to_order(self):
        source = _atom_array(
            chain_ids=["B", "A"],
            res_ids=[1, 1],
            atom_names=["N", "N"],
            res_names=["GLY", "ALA"],
            coords=[[2, 0, 0], [1, 0, 0]],
        )
        internal = _atom_array(
            chain_ids=["A", "B"],
            res_ids=[1, 1],
            atom_names=["N", "N"],
            res_names=["ALA", "GLY"],
            coords=[[1, 0, 0], [2, 0, 0]],
        )

        chain_id_map = build_chain_id_map(source, internal, prefer_identity=False)

        self.assertEqual(chain_id_map["internal_to_source"], {"A": "B", "B": "A"})

    def test_map_source_coords_to_internal_with_reference_fallback(self):
        source = _atom_array(
            chain_ids=["A", "A", "B"],
            res_ids=[1, 1, 1],
            atom_names=["N", "CA", "C1"],
            res_names=["ALA", "ALA", "ATP"],
            coords=[[1, 0, 0], [2, 0, 0], [3, 0, 0]],
        )
        internal = _atom_array(
            chain_ids=["X", "X", "Y", "Y"],
            res_ids=[1, 1, 1, 1],
            atom_names=["N", "CA", "C1", "O1"],
            res_names=["ALA", "ALA", "ATP", "ATP"],
            coords=[[9, 0, 0], [9, 1, 0], [9, 2, 0], [9, 3, 0]],
        )
        chain_id_map = build_chain_id_map(source, internal)
        coords, missing_atoms = map_source_coords_to_internal(
            source_coord_map=build_source_coord_map(source),
            internal_atom_array=internal,
            chain_id_map=chain_id_map,
            missing_atom_policy="reference",
        )

        self.assertEqual(len(missing_atoms), 1)
        self.assertTrue(np.allclose(coords[0], [1, 0, 0]))
        self.assertTrue(np.allclose(coords[1], [2, 0, 0]))
        self.assertTrue(np.allclose(coords[2], [3, 0, 0]))
        self.assertTrue(np.allclose(coords[3], [9, 3, 0]))

    def test_map_source_coords_to_internal_errors_on_missing_atoms(self):
        source = _atom_array(
            chain_ids=["A"],
            res_ids=[1],
            atom_names=["N"],
            res_names=["ALA"],
            coords=[[1, 0, 0]],
        )
        internal = _atom_array(
            chain_ids=["A", "A"],
            res_ids=[1, 1],
            atom_names=["N", "CA"],
            res_names=["ALA", "ALA"],
            coords=[[1, 0, 0], [2, 0, 0]],
        )
        chain_id_map = build_chain_id_map(source, internal)

        with self.assertRaises(RuntimeError):
            map_source_coords_to_internal(
                source_coord_map=build_source_coord_map(source),
                internal_atom_array=internal,
                chain_id_map=chain_id_map,
                missing_atom_policy="error",
            )


if __name__ == "__main__":
    unittest.main()
