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

import numpy as np
import pytest
from biotite.structure import AtomArray, BondList, BondType

from protenix.data.core.oxygen_branches import absent_oxygen_branch


def fixture():
    comp = AtomArray(7)
    comp.atom_name = np.array(["O5'", "P", "OP1", "OP2", "OP3", "C5'", "C4'"])
    comp.element = np.array(["O", "P", "O", "O", "O", "C", "C"])
    comp.res_name[:] = "G"
    comp.set_annotation("charge", np.zeros(7, dtype=int))
    comp.bonds = BondList(
        7, np.array([[0, 1, 1], [1, 2, 2], [1, 3, 1], [1, 4, 1], [0, 5, 1], [5, 6, 1]])
    )
    raw = comp[[0, 5, 6]].copy()
    raw.chain_id[:] = "A"
    raw.res_id[:] = 1
    raw.set_annotation("label_entity_id", np.array(["1"] * 3))
    ext = raw[:1].copy()
    ext.res_id[:] = 2
    ext.res_name[:] = "CAP"
    ext.atom_name[:] = "PG"
    ext.element[:] = "P"
    raw += ext
    raw.bonds = BondList(4, np.array([[0, 1, 1], [1, 2, 1], [0, 3, 1]]))
    return comp, raw


def test_missing_phosphate_branch():
    comp, raw = fixture()
    assert absent_oxygen_branch(raw, comp, "O5'", [0]) == ["OP1", "OP2", "OP3", "P"]


@pytest.mark.parametrize(
    "negative",
    [
        "observed",
        "other_copy",
        "missing_center_copy",
        "double",
        "metal",
        "unknown_order",
        "coordination",
        "center_ring",
        "multiple",
        "ring",
        "unanchored",
        "charged",
        "ccd_double",
    ],
)
def test_ambiguous_or_protected_branch_retained(negative):
    comp, raw = fixture()
    indices = [0]
    if negative in {"observed", "other_copy", "missing_center_copy"}:
        extra = raw[:1].copy()
        extra.atom_name[:] = "OP1"
        extra.element[:] = "O"
        if negative != "observed":
            extra.chain_id[:] = "B"
        raw += extra
        if negative == "other_copy":
            copy = raw[:4].copy()
            copy.chain_id[:] = "B"
            raw += copy
            indices.append(5)
    elif negative == "double":
        raw.bonds = BondList(4, np.array([[0, 1, 1], [1, 2, 1], [0, 3, 2]]))
    elif negative in {"unknown_order", "coordination"}:
        order = BondType.ANY if negative == "unknown_order" else BondType.COORDINATION
        raw.bonds = BondList(4, np.array([[0, 1, 1], [1, 2, 1], [0, 3, order]]))
    elif negative == "center_ring":
        comp.bonds.add_bond(1, 5, BondType.SINGLE)
    elif negative == "metal":
        raw.element[3] = "FE"
    elif negative == "multiple":
        raw += raw[3:4].copy()
        raw.bonds = BondList(5, np.array([[0, 1, 1], [1, 2, 1], [0, 3, 1], [0, 4, 1]]))
    elif negative == "ring":
        comp.bonds.add_bond(2, 3, BondType.SINGLE)
    elif negative == "unanchored":
        raw.atom_name[1:3] = ["REMOTE1", "REMOTE2"]
    elif negative == "charged":
        comp.charge[0] = 1
    elif negative == "ccd_double":
        comp.bonds.add_bond(0, 1, BondType.DOUBLE)
    assert absent_oxygen_branch(raw, comp, "O5'", indices) == []


def test_disconnected_observation_cannot_anchor_branch():
    comp, raw = fixture()
    comp.bonds = BondList(
        7, np.array([[0, 1, 1], [1, 2, 2], [1, 3, 1], [1, 4, 1], [0, 5, 1]])
    )
    raw.atom_name[1] = "REMOTE"
    assert absent_oxygen_branch(raw, comp, "O5'", [0]) == []


def test_named_confidence_slots_survive_leading_removal():
    from protenix.data.core.parser import AddAtomArrayAnnot

    comp, _ = fixture()
    comp.set_annotation("mol_type", np.array(["rna"] * len(comp)))
    original = AddAtomArrayAnnot.add_tokatom_idx(comp.copy())
    kept = np.array([0, 5, 6])
    pruned = AddAtomArrayAnnot.add_tokatom_idx(comp[kept])
    assert np.array_equal(pruned.tokatom_idx, original.tokatom_idx[kept])
