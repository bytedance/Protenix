import unittest

from protenix.utils.epitope_alias import apply_epitopes_alias


class TestEpitopeAlias(unittest.TestCase):
    def test_epitopes_positions_alias_to_pocket(self):
        sample = {
            "name": "x",
            "sequences": [],
            "epitopes": {
                "binder_chain": {"entity": 2, "copy": 1},
                "antigen": {"entity": 1, "copy": 1},
                "positions": [45, 46],
                "max_distance": 7.5,
            },
        }
        out = apply_epitopes_alias(sample)

        self.assertNotIn("epitopes", out)
        self.assertIn("constraint", out)
        self.assertIn("pocket", out["constraint"])

        pocket = out["constraint"]["pocket"]
        self.assertEqual(pocket["binder_chain"], {"entity": 2, "copy": 1})
        self.assertEqual(pocket["max_distance"], 7.5)
        self.assertEqual(
            pocket["contact_residues"],
            [
                {"entity": 1, "copy": 1, "position": 45},
                {"entity": 1, "copy": 1, "position": 46},
            ],
        )

    def test_explicit_pocket_wins_over_epitopes(self):
        sample = {
            "name": "x",
            "sequences": [],
            "constraint": {
                "pocket": {
                    "binder_chain": {"entity": 2, "copy": 1},
                    "contact_residues": [{"entity": 1, "copy": 1, "position": 69}],
                    "max_distance": 8,
                }
            },
            "epitopes": {
                "binder_chain": {"entity": 2, "copy": 1},
                "antigen": {"entity": 1, "copy": 1},
                "positions": [45, 46],
                "max_distance": 6,
            },
        }
        out = apply_epitopes_alias(sample)

        self.assertNotIn("epitopes", out)
        self.assertEqual(out["constraint"]["pocket"]["contact_residues"][0]["position"], 69)

