"""Tests for MSA encoding vectorized sequence parsing.

Verifies the numpy LUT-based approach produces identical msa_arr and del_arr
as the original per-character Python loop for protein/DNA/RNA sequences.
"""

import unittest

import numpy as np


# ---------- Reference (upstream) implementation ----------
def _encode_msa_old(sequences, cmap):
    """Original per-character loop implementation."""
    if not sequences:
        return np.zeros((0, 0), dtype=np.int32), np.zeros((0, 0), dtype=np.int32)

    rows, cols = len(sequences), sum(1 for c in sequences[0] if c in cmap)
    msa_arr = np.zeros((rows, cols), dtype=np.int32)
    del_arr = np.zeros((rows, cols), dtype=np.int32)

    for i, seq in enumerate(sequences):
        d_count, c_idx = 0, 0
        for char in seq:
            val = cmap.get(char, -1)
            if val == -1:
                d_count += 1
            else:
                if c_idx < cols:
                    msa_arr[i, c_idx] = val
                    del_arr[i, c_idx] = d_count
                d_count, c_idx = 0, c_idx + 1
    return msa_arr, del_arr


# ---------- New (vectorized) implementation ----------
def _encode_msa_new(sequences, cmap):
    """Numpy LUT-based implementation."""
    if not sequences:
        return np.zeros((0, 0), dtype=np.int32), np.zeros((0, 0), dtype=np.int32)

    lut = np.full(128, -1, dtype=np.int32)
    for char, val in cmap.items():
        lut[ord(char)] = val

    cols = sum(1 for c in sequences[0] if c in cmap)
    rows = len(sequences)
    msa_arr = np.zeros((rows, cols), dtype=np.int32)
    del_arr = np.zeros((rows, cols), dtype=np.int32)

    for i, seq in enumerate(sequences):
        char_arr = np.frombuffer(seq.encode("ascii"), dtype=np.uint8)
        mapped = lut[char_arr]
        aligned_mask = mapped >= 0
        aligned_vals = mapped[aligned_mask]
        cum_inserts = np.cumsum(~aligned_mask)
        insert_before_aligned = cum_inserts[aligned_mask]
        del_counts = np.diff(insert_before_aligned, prepend=0)
        n = min(len(aligned_vals), cols)
        msa_arr[i, :n] = aligned_vals[:n]
        del_arr[i, :n] = del_counts[:n]
    return msa_arr, del_arr


# Protein char map (from protenix)
PROTEIN_CMAP = {
    'A': 0, 'B': 3, 'C': 4, 'D': 3, 'E': 6, 'F': 13, 'G': 7, 'H': 8,
    'I': 9, 'J': 20, 'K': 11, 'L': 10, 'M': 12, 'N': 2, 'O': 20, 'P': 14,
    'Q': 5, 'R': 1, 'S': 15, 'T': 16, 'U': 4, 'V': 19, 'W': 17, 'X': 20,
    'Y': 18, 'Z': 6, '-': 31,
}


class TestMSAEncodingEquivalence(unittest.TestCase):
    """Compare old (loop) vs new (LUT) MSA encoding."""

    def _assert_equal(self, sequences, cmap):
        msa_old, del_old = _encode_msa_old(sequences, cmap)
        msa_new, del_new = _encode_msa_new(sequences, cmap)
        np.testing.assert_array_equal(msa_old, msa_new, err_msg="msa_arr mismatch")
        np.testing.assert_array_equal(del_old, del_new, err_msg="del_arr mismatch")

    def test_simple_aligned(self):
        """No insertions, pure aligned positions."""
        seqs = ["ACDEFG", "GHIKLM"]
        self._assert_equal(seqs, PROTEIN_CMAP)

    def test_with_gaps(self):
        """Gaps (-) are aligned positions with value 31."""
        seqs = ["AC-EFG", "GH-KLM"]
        self._assert_equal(seqs, PROTEIN_CMAP)

    def test_with_insertions(self):
        """Lowercase chars are insertions (not in cmap), should count as deletions."""
        # In real MSAs, insertions are lowercase letters not in the cmap
        # We simulate by adding chars that don't exist in cmap
        cmap = {'A': 0, 'C': 1, 'G': 2, 'T': 3, '-': 4}
        seqs = [
            "AxxCG",   # 2 insertions between A and C
            "ACxxG",   # 2 insertions between C and G
        ]
        self._assert_equal(seqs, cmap)

    def test_multiple_insertion_runs(self):
        """Multiple runs of insertions in one sequence."""
        cmap = {'A': 0, 'C': 1, 'G': 2, '-': 3}
        seqs = [
            "AxCxxG",   # 1 insert before C, 2 before G
            "AxxCxG",   # 2 before C, 1 before G
        ]
        self._assert_equal(seqs, cmap)

    def test_leading_insertions(self):
        """Insertions at the very start of sequence."""
        cmap = {'A': 0, 'C': 1, 'G': 2, '-': 3}
        seqs = [
            "xxACG",   # 2 leading insertions
            "ACG",
        ]
        self._assert_equal(seqs, cmap)

    def test_trailing_insertions(self):
        """Insertions at the end (should be ignored, beyond aligned columns)."""
        cmap = {'A': 0, 'C': 1, 'G': 2, '-': 3}
        seqs = [
            "ACGxx",
            "ACG",
        ]
        self._assert_equal(seqs, cmap)

    def test_empty_sequences(self):
        self._assert_equal([], PROTEIN_CMAP)

    def test_single_column(self):
        seqs = ["A", "G", "C"]
        self._assert_equal(seqs, PROTEIN_CMAP)

    def test_all_gaps(self):
        seqs = ["---", "---"]
        self._assert_equal(seqs, PROTEIN_CMAP)

    def test_real_protein_sequences(self):
        """Longer realistic protein MSA with insertions."""
        seqs = [
            "MAQGSHQIDFQVLHD",
            "MAQGxSHQIDxxFQVLHD",
            "MAQGSHxQIDFQVxLHD",
        ]
        # First seq defines cols (all uppercase = 16 columns)
        self._assert_equal(seqs, PROTEIN_CMAP)

    def test_full_protein_cmap(self):
        """Exercise all protein chars."""
        seqs = ["ABCDEFGHIJKLMNOPQRSTUVWXYZ-"]
        self._assert_equal(seqs, PROTEIN_CMAP)

    def test_large_msa(self):
        """Stress test: 1000 sequences of length 200."""
        np.random.seed(42)
        chars = list(PROTEIN_CMAP.keys()) + list("abcdefgh")  # mix aligned + insertion
        seqs = []
        # First sequence defines columns - use only aligned chars
        aligned_chars = list(PROTEIN_CMAP.keys())
        seqs.append("".join(np.random.choice(aligned_chars, 200)))
        # Rest can have insertions
        for _ in range(999):
            seqs.append("".join(np.random.choice(chars, 200)))
        self._assert_equal(seqs, PROTEIN_CMAP)


class TestMSAEncodingIntegration(unittest.TestCase):
    """Test via the actual MSACore.sequences_to_array method."""

    def test_sequences_to_array_protein(self):
        from protenix.data.msa.msa_utils import MSACore, PROTEIN_CHAIN
        seqs = ["ACDEFGHIKLMNPQRSTVWY", "ACDEFGHIKLMNPQRSTVWY"]
        msa, dele = MSACore.sequences_to_array(seqs, chain_type=PROTEIN_CHAIN)
        self.assertEqual(msa.shape, (2, 20))
        self.assertEqual(dele.shape, (2, 20))
        # All positions aligned, no deletions
        np.testing.assert_array_equal(dele, np.zeros_like(dele))

    def test_sequences_to_array_empty(self):
        from protenix.data.msa.msa_utils import MSACore, PROTEIN_CHAIN
        msa, dele = MSACore.sequences_to_array([], chain_type=PROTEIN_CHAIN)
        self.assertEqual(msa.shape, (0, 0))


if __name__ == "__main__":
    unittest.main()
