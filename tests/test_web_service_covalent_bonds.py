import unittest

from protenix.web_service.viewer import ProtenixInputViewer


class TestProtenixInputViewerCovalentBonds(unittest.TestCase):
    def test_polymer_position_uses_sequence_length(self):
        viewer = ProtenixInputViewer()
        viewer.name.value = "demo"

        viewer.add_dna_rna_protein_callback(None)
        protein = viewer.dna_rna_protein_entities[0].children[0]
        protein.sequence_text.value = "AAAAA"

        viewer.add_ligand_smiles_callback(None)
        ligand = viewer.ligand_smiles_entities[0].children[0]
        ligand.smiles_text.value = "C"

        viewer.covalent_bonds.on_add_convalent_bond(None)
        bond = viewer.covalent_bonds.convalent_bonds[0]
        bond.children[1].value = 5
        bond.children[2].value = "CA"
        bond.children[3].value = 2
        bond.children[4].value = 1
        bond.children[5].value = "C1"

        self.assertEqual(
            viewer.get_result()["covalent_bonds"],
            [
                {
                    "left_entity": 1,
                    "left_position": 5,
                    "left_atom": "CA",
                    "right_entity": 2,
                    "right_position": 1,
                    "right_atom": "C1",
                }
            ],
        )

        bond.children[1].value = 6
        with self.assertRaisesRegex(AssertionError, "left position out of range"):
            viewer.get_result()


if __name__ == "__main__":
    unittest.main()
