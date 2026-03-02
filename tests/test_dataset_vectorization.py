"""Tests for Dataset vectorize pandas/numpy operations.

Verifies that vectorized pandas/numpy ops produce identical results to
the original apply/lambda/vectorize patterns.
"""

import unittest

import numpy as np
import pandas as pd


class TestMolGroupType(unittest.TestCase):
    """Verify vectorized mol_group_type matches df.apply(sorted join)."""

    def _old_mol_group_type(self, df):
        return df.apply(
            lambda row: "_".join(
                sorted([str(row["mol_1_type"]), str(row["mol_2_type"]).replace("nan", "intra")])
            ),
            axis=1,
        )

    def _new_mol_group_type(self, df):
        col1 = df["mol_1_type"].astype(str)
        col2 = df["mol_2_type"].astype(str).str.replace("nan", "intra", regex=False)
        lo = np.where(col1.values <= col2.values, col1.values, col2.values)
        hi = np.where(col1.values <= col2.values, col2.values, col1.values)
        return pd.Series(np.char.add(np.char.add(lo, "_"), hi), index=df.index)

    def test_basic(self):
        df = pd.DataFrame({
            "mol_1_type": ["protein", "dna", "ligand", "protein"],
            "mol_2_type": ["ligand", "protein", "protein", float("nan")],
        })
        old = self._old_mol_group_type(df)
        new = self._new_mol_group_type(df)
        pd.testing.assert_series_equal(old, new, check_names=False)

    def test_all_nan(self):
        df = pd.DataFrame({
            "mol_1_type": ["protein", "dna"],
            "mol_2_type": [float("nan"), float("nan")],
        })
        old = self._old_mol_group_type(df)
        new = self._new_mol_group_type(df)
        pd.testing.assert_series_equal(old, new, check_names=False)

    def test_same_types(self):
        df = pd.DataFrame({
            "mol_1_type": ["protein", "protein"],
            "mol_2_type": ["protein", "protein"],
        })
        old = self._old_mol_group_type(df)
        new = self._new_mol_group_type(df)
        pd.testing.assert_series_equal(old, new, check_names=False)

    def test_sorting_order(self):
        """Verify sorted order is consistent."""
        df = pd.DataFrame({
            "mol_1_type": ["zzz", "aaa"],
            "mol_2_type": ["aaa", "zzz"],
        })
        old = self._old_mol_group_type(df)
        new = self._new_mol_group_type(df)
        pd.testing.assert_series_equal(old, new, check_names=False)
        # Both should be "aaa_zzz"
        self.assertTrue(all(v == "aaa_zzz" for v in new))


class TestGetContiguousArray(unittest.TestCase):
    """Verify np.searchsorted matches dict-based vectorize for contiguous remapping."""

    def _old_get_contiguous(self, array):
        array_uniq = np.sort(np.unique(array))
        map_dict = {i: idx for idx, i in enumerate(array_uniq)}
        return np.vectorize(map_dict.get)(array)

    def _new_get_contiguous(self, array):
        array_uniq = np.sort(np.unique(array))
        return np.searchsorted(array_uniq, array)

    def test_basic(self):
        arr = np.array([0, 1, 2, 18, 20])
        np.testing.assert_array_equal(
            self._old_get_contiguous(arr),
            self._new_get_contiguous(arr),
        )

    def test_already_contiguous(self):
        arr = np.array([0, 1, 2, 3, 4])
        np.testing.assert_array_equal(
            self._old_get_contiguous(arr),
            self._new_get_contiguous(arr),
        )

    def test_with_gaps(self):
        arr = np.array([0, 0, 0, 5, 5, 10, 10, 10, 20])
        np.testing.assert_array_equal(
            self._old_get_contiguous(arr),
            self._new_get_contiguous(arr),
        )

    def test_reversed(self):
        arr = np.array([20, 18, 2, 1, 0])
        np.testing.assert_array_equal(
            self._old_get_contiguous(arr),
            self._new_get_contiguous(arr),
        )

    def test_single_value(self):
        arr = np.array([42, 42, 42])
        np.testing.assert_array_equal(
            self._old_get_contiguous(arr),
            self._new_get_contiguous(arr),
        )


class TestShuffleViaSearchsorted(unittest.TestCase):
    """Verify searchsorted-based shuffle matches dict-based vectorize."""

    def _old_shuffle(self, x, seed):
        np.random.seed(seed)
        x_unique = np.sort(np.unique(x))
        x_shuffled = x_unique.copy()
        np.random.shuffle(x_shuffled)
        map_dict = dict(zip(x_unique, x_shuffled))
        return np.vectorize(map_dict.get)(x).copy()

    def _new_shuffle(self, x, seed):
        np.random.seed(seed)
        x_unique = np.sort(np.unique(x))
        x_shuffled = x_unique.copy()
        np.random.shuffle(x_shuffled)
        indices = np.searchsorted(x_unique, x)
        return x_shuffled[indices].copy()

    def test_basic(self):
        arr = np.array([0, 0, 1, 1, 2, 2, 3])
        for seed in [42, 123, 999]:
            np.testing.assert_array_equal(
                self._old_shuffle(arr, seed),
                self._new_shuffle(arr, seed),
            )

    def test_single_value(self):
        arr = np.array([5, 5, 5])
        np.testing.assert_array_equal(
            self._old_shuffle(arr, 42),
            self._new_shuffle(arr, 42),
        )

    def test_large_array(self):
        np.random.seed(0)
        arr = np.random.randint(0, 50, size=10000)
        np.testing.assert_array_equal(
            self._old_shuffle(arr, 42),
            self._new_shuffle(arr, 42),
        )


class TestIsinVsApplyLambda(unittest.TestCase):
    """Verify df.isin matches df.apply(lambda x: x in set)."""

    def test_basic(self):
        valid_set = {"chain", "interface", "ligand"}
        df = pd.DataFrame({
            "eval_type": ["chain", "other", "interface", "unknown", "ligand"]
        })
        old = df["eval_type"].apply(lambda x: x in valid_set)
        new = df["eval_type"].isin(valid_set)
        pd.testing.assert_series_equal(old, new)

    def test_empty_df(self):
        valid_set = {"chain"}
        df = pd.DataFrame({"eval_type": pd.Series([], dtype=str)})
        old = df["eval_type"].apply(lambda x: x in valid_set)
        new = df["eval_type"].isin(valid_set)
        # Both are empty; apply returns object dtype on empty, isin returns bool
        # The actual filtering result (empty mask) is identical
        self.assertEqual(len(old), 0)
        self.assertEqual(len(new), 0)


class TestCalcWeightsSortedEntityId(unittest.TestCase):
    """Verify vectorized pdb_sorted_entity_id matches apply-based version."""

    def _old_sorted_entity(self, df):
        return df.apply(
            lambda x: f"{x['pdb_id']}_{x['assembly_id']}_{'_'.join(sorted([str(x['entity_1_id']), str(x['entity_2_id'])]))}",
            axis=1,
        )

    def _new_sorted_entity(self, df):
        e1 = df["entity_1_id"].astype(str).values
        e2 = df["entity_2_id"].astype(str).values
        lo = np.where(e1 <= e2, e1, e2)
        hi = np.where(e1 <= e2, e2, e1)
        return (
            df["pdb_id"].astype(str) + "_"
            + df["assembly_id"].astype(str) + "_"
            + lo + "_" + hi
        )

    def test_basic(self):
        df = pd.DataFrame({
            "pdb_id": ["1abc", "1abc", "2def"],
            "assembly_id": [1, 1, 2],
            "entity_1_id": [1, 3, 2],
            "entity_2_id": [3, 1, 1],
        })
        old = self._old_sorted_entity(df)
        new = self._new_sorted_entity(df)
        pd.testing.assert_series_equal(old, new, check_names=False)

    def test_with_nan(self):
        df = pd.DataFrame({
            "pdb_id": ["1abc"],
            "assembly_id": [1],
            "entity_1_id": [1],
            "entity_2_id": [float("nan")],
        })
        old = self._old_sorted_entity(df)
        new = self._new_sorted_entity(df)
        pd.testing.assert_series_equal(old, new, check_names=False)

    def test_map_vs_apply_for_member_count(self):
        """Verify .map(dict) matches .apply(lambda: dict[key])."""
        df = pd.DataFrame({
            "key": ["a", "b", "a", "c", "b", "a"],
        })
        count_dict = {"a": 3, "b": 2, "c": 1}
        old = df.apply(lambda x: count_dict[x["key"]], axis=1)
        new = df["key"].map(count_dict)
        pd.testing.assert_series_equal(old, new, check_names=False)


if __name__ == "__main__":
    unittest.main()
