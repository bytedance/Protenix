"""Tests for PR 1: featurizer.py optimizations.

Verifies that:
1. torch.from_numpy produces same dtype/shape/values as torch.Tensor().long()
2. ndarray.copy() matches copy.deepcopy for simple arrays
3. Pre-allocated ref_pos matches append+concatenate
4. _get_atom_to_token_idx caching + vectorized bond adjacency
5. Vectorized near_atoms (np.isin) matches list comprehension
"""

import copy
import unittest

import numpy as np
import torch


class TestTorchFromNumpyEquivalence(unittest.TestCase):
    """Verify torch.from_numpy(arr.astype(np.int64)) == torch.Tensor(arr).long()"""

    def _check_equivalence(self, arr):
        old = torch.Tensor(arr).long()
        new = torch.from_numpy(np.asarray(arr).astype(np.int64))
        self.assertEqual(old.dtype, new.dtype, f"dtype mismatch: {old.dtype} vs {new.dtype}")
        self.assertEqual(old.shape, new.shape, f"shape mismatch: {old.shape} vs {new.shape}")
        self.assertTrue(torch.equal(old, new), "values differ")

    def test_1d_int_array(self):
        self._check_equivalence(np.array([0, 1, 5, 100, -3]))

    def test_1d_float_array_truncates(self):
        """float values get truncated to int64 the same way"""
        arr = np.array([0.0, 1.9, 2.1, -0.5])
        old = torch.Tensor(arr).long()
        new = torch.from_numpy(arr.astype(np.int64))
        self.assertTrue(torch.equal(old, new))

    def test_2d_array(self):
        self._check_equivalence(np.array([[1, 2], [3, 4], [5, 6]]))

    def test_bool_array(self):
        """is_resolved, is_protein, etc. are bool arrays"""
        arr = np.array([True, False, True, True, False])
        self._check_equivalence(arr.astype(int))

    def test_empty_array(self):
        self._check_equivalence(np.array([], dtype=int))

    def test_large_values_precision_improvement(self):
        """torch.Tensor(arr).long() goes via float32, losing precision for large ints.
        torch.from_numpy preserves exact values. This is a correctness improvement.
        In practice, featurizer arrays contain small IDs so this doesn't affect results,
        but from_numpy is strictly more correct."""
        arr = np.array([0, 2**30, 2**31 - 1, -2**31])
        old = torch.Tensor(arr).long()
        new = torch.from_numpy(arr.astype(np.int64))
        # 2**31 - 1 rounds to 2**31 via float32
        self.assertNotEqual(old[2].item(), 2**31 - 1)  # old is WRONG
        self.assertEqual(new[2].item(), 2**31 - 1)      # new is CORRECT

    def test_typical_featurizer_values(self):
        """Values actually seen in featurizer: small IDs, masks, indices."""
        for arr in [
            np.arange(500),                          # token indices
            np.array([0, 0, 0, 1, 1, 2, 2, 2, 3]),  # asym_id
            np.array([1, 0, 1, 1, 0]),               # masks
        ]:
            self._check_equivalence(arr)


class TestNdarrayCopyVsDeepCopy(unittest.TestCase):
    """Verify ndarray.copy() == copy.deepcopy(ndarray) for string arrays."""

    def test_string_array_copy(self):
        arr = np.array(["CA", "N", "C", "O", "CB"])
        old = copy.deepcopy(arr)
        new = arr.copy()
        np.testing.assert_array_equal(old, new)

    def test_mutation_independence(self):
        arr = np.array(["CA", "N", "C"])
        new = arr.copy()
        new[0] = "XX"
        self.assertEqual(arr[0], "CA")  # original unchanged


class TestRefPosPreallocation(unittest.TestCase):
    """Verify pre-allocated ref_pos matches append+concatenate."""

    def test_mask_assign_vs_append_concat(self):
        np.random.seed(42)
        N = 50
        ref_pos = np.random.randn(N, 3).astype(np.float32)
        ref_space_uid = np.array([0] * 20 + [1] * 15 + [2] * 15)

        # Old: append + concat
        result_old = []
        for uid in np.unique(ref_space_uid):
            chunk = ref_pos[ref_space_uid == uid].copy()
            chunk -= chunk.mean(axis=0)  # simplified "random_transform" with centralize
            result_old.append(chunk)
        result_old = np.concatenate(result_old)

        # New: pre-allocate + mask assign
        result_new = np.empty_like(ref_pos)
        for uid in np.unique(ref_space_uid):
            mask = ref_space_uid == uid
            chunk = ref_pos[mask].copy()
            chunk -= chunk.mean(axis=0)
            result_new[mask] = chunk

        np.testing.assert_allclose(result_old, result_new, atol=1e-6)


class TestAtomToTokenIdx(unittest.TestCase):
    """Test _get_atom_to_token_idx logic and caching."""

    def _build_mapping_old(self, tokens, n_atoms):
        """Original: dict comprehension → list comprehension"""
        atom_to_token_idx_dict = {}
        for idx, token in enumerate(tokens):
            for atom_idx in token["atom_indices"]:
                atom_to_token_idx_dict[atom_idx] = idx
        return [atom_to_token_idx_dict[i] for i in range(n_atoms)]

    def _build_mapping_new(self, tokens, n_atoms):
        """New: numpy array with -1 init"""
        atom_idx_to_token_idx = np.full(n_atoms, -1, dtype=int)
        for idx, token in enumerate(tokens):
            for atom_idx in token["atom_indices"]:
                atom_idx_to_token_idx[atom_idx] = idx
        return atom_idx_to_token_idx

    def test_equivalence(self):
        tokens = [
            {"atom_indices": [0, 1, 2, 3]},
            {"atom_indices": [4, 5]},
            {"atom_indices": [6, 7, 8]},
            {"atom_indices": [9]},
        ]
        n_atoms = 10
        old = self._build_mapping_old(tokens, n_atoms)
        new = self._build_mapping_new(tokens, n_atoms)
        np.testing.assert_array_equal(old, new)

    def test_single_atom_tokens(self):
        tokens = [{"atom_indices": [i]} for i in range(5)]
        old = self._build_mapping_old(tokens, 5)
        new = self._build_mapping_new(tokens, 5)
        np.testing.assert_array_equal(old, new)


class TestVectorizedBondAdjacency(unittest.TestCase):
    """Verify vectorized indexing matches the Python loop for bond adjacency."""

    def test_equivalence(self):
        np.random.seed(42)
        num_tokens = 8
        atom_idx_to_token_idx = np.array([0, 0, 0, 1, 1, 2, 2, 3, 4, 4, 5, 6, 6, 7])
        # Random bonds
        kept_bonds = np.array([
            [0, 5], [1, 6], [3, 8], [9, 12], [2, 13], [4, 11],
        ])

        # Old: Python loop
        old_matrix = np.zeros((num_tokens, num_tokens), dtype=int)
        bond_token_i = atom_idx_to_token_idx[kept_bonds[:, 0]]
        bond_token_j = atom_idx_to_token_idx[kept_bonds[:, 1]]
        for i, j in zip(bond_token_i, bond_token_j):
            old_matrix[i, j] = 1
            old_matrix[j, i] = 1

        # New: vectorized indexing
        new_matrix = np.zeros((num_tokens, num_tokens), dtype=int)
        bt_i = atom_idx_to_token_idx[kept_bonds[:, 0]]
        bt_j = atom_idx_to_token_idx[kept_bonds[:, 1]]
        new_matrix[bt_i, bt_j] = 1
        new_matrix[bt_j, bt_i] = 1

        np.testing.assert_array_equal(old_matrix, new_matrix)

    def test_duplicate_bonds_same_token_pair(self):
        """Multiple bonds between same token pair should still be 1."""
        num_tokens = 4
        atom_idx_to_token_idx = np.array([0, 0, 1, 1, 2, 3])
        kept_bonds = np.array([[0, 2], [1, 3]])  # both map to token (0,1)

        old_matrix = np.zeros((num_tokens, num_tokens), dtype=int)
        for bond in kept_bonds:
            i, j = atom_idx_to_token_idx[bond[0]], atom_idx_to_token_idx[bond[1]]
            old_matrix[i, j] = 1
            old_matrix[j, i] = 1

        new_matrix = np.zeros((num_tokens, num_tokens), dtype=int)
        bt_i = atom_idx_to_token_idx[kept_bonds[:, 0]]
        bt_j = atom_idx_to_token_idx[kept_bonds[:, 1]]
        new_matrix[bt_i, bt_j] = 1
        new_matrix[bt_j, bt_i] = 1

        np.testing.assert_array_equal(old_matrix, new_matrix)


class TestVectorizedNearAtoms(unittest.TestCase):
    """Verify np.isin vectorized approach matches list comprehension."""

    def test_equivalence(self):
        n = 20
        near_atom_indices = np.array([2, 5, 7, 10, 15])
        is_resolved = np.random.randint(0, 2, size=n).astype(bool)

        # Old: list comprehension
        old = [
            True if ((i in near_atom_indices) and is_resolved[i]) else False
            for i in range(n)
        ]

        # New: vectorized
        new = (
            np.isin(np.arange(n), near_atom_indices)
            & is_resolved.astype(bool)
        )

        np.testing.assert_array_equal(old, new)

    def test_empty_indices(self):
        n = 10
        near_atom_indices = np.array([], dtype=int)
        is_resolved = np.ones(n, dtype=bool)

        old = [
            True if ((i in near_atom_indices) and is_resolved[i]) else False
            for i in range(n)
        ]
        new = np.isin(np.arange(n), near_atom_indices) & is_resolved
        np.testing.assert_array_equal(old, new)

    def test_all_resolved_all_near(self):
        n = 5
        near_atom_indices = np.arange(n)
        is_resolved = np.ones(n, dtype=bool)

        old = [True for _ in range(n)]
        new = np.isin(np.arange(n), near_atom_indices) & is_resolved
        np.testing.assert_array_equal(old, new)


if __name__ == "__main__":
    unittest.main()
