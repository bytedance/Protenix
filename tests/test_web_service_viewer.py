import unittest

from protenix.web_service.viewer import DnaRnaProteinEntityWidget


class TestDnaRnaProteinEntityWidget(unittest.TestCase):
    def _get_modification(self, molecule_type):
        widget = DnaRnaProteinEntityWidget()
        widget.molecule_type_dropdown.value = molecule_type
        widget.sequence_text.value = "A"
        widget.add_modification_callback(None)

        entity_key = {
            "Protein": "proteinChain",
            "DNA": "dnaSequence",
            "RNA": "rnaSequence",
        }[molecule_type]
        return widget.get_result()[entity_key]["modifications"][0]

    def test_protein_modification_uses_inference_schema_keys(self):
        self.assertEqual(
            self._get_modification("Protein"),
            {"ptmType": "CCD_XXX", "ptmPosition": 1},
        )

    def test_dna_modification_uses_inference_schema_keys(self):
        self.assertEqual(
            self._get_modification("DNA"),
            {"modificationType": "CCD_XXX", "basePosition": 1},
        )

    def test_rna_modification_uses_inference_schema_keys(self):
        self.assertEqual(
            self._get_modification("RNA"),
            {"modificationType": "CCD_XXX", "basePosition": 1},
        )


if __name__ == "__main__":
    unittest.main()
