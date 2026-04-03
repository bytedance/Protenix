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


from collections.abc import Sequence

import numpy as np
from biotite.file import InvalidFileError
from biotite.structure import AtomArray, AtomArrayStack
from biotite.structure.bonds import BondList, BondType
from biotite.structure.box import coord_to_fraction, fraction_to_coord
from biotite.structure.io.pdbx import convert as pdbx_convert

PDBX_BOND_TYPE_ID_TO_TYPE = {
    # Although a covalent bond, could in theory have a higher bond order,
    # practically inter-residue bonds are always single
    "covale": BondType.SINGLE,
    "covale_base": BondType.SINGLE,
    "covale_phosphate": BondType.SINGLE,
    "covale_sugar": BondType.SINGLE,
    "disulf": BondType.SINGLE,
    "modres": BondType.SINGLE,
    "modres_link": BondType.SINGLE,
    # "metalc": BondType.COORDINATION,  # no metalc for Protenix
}


PDBX_BOND_TYPE_TO_ORDER = {
    BondType.SINGLE: "sing",
    BondType.DOUBLE: "doub",
    BondType.TRIPLE: "trip",
    BondType.QUADRUPLE: "quad",
    BondType.AROMATIC_SINGLE: "sing",
    BondType.AROMATIC_DOUBLE: "doub",
    BondType.AROMATIC_TRIPLE: "trip",
    # These are masked later, it is merely added here to avoid a KeyError
    BondType.ANY: "",
    # BondType.COORDINATION: "",
}


FIND_MATCHES_SWITCH_THRES = 4000000


def _find_matches_by_dense_array(query_arrays, reference_arrays):
    """
    For each index in the `query_arrays` find the indices in the
    `reference_arrays` where all query values match the reference counterpart.
    If no match is found for a query, the corresponding index is -1.
    """
    match_masks_for_all_columns = np.stack(
        [
            query[:, np.newaxis] == reference[np.newaxis, :]
            for query, reference in zip(query_arrays, reference_arrays)
        ],
        axis=-1,
    )
    match_masks = np.all(match_masks_for_all_columns, axis=-1)
    query_matches, reference_matches = np.where(match_masks)

    # Duplicate matches indicate that an atom from the query cannot
    # be uniquely matched to an atom in the reference
    unique_query_matches, counts = np.unique(query_matches, return_counts=True)
    if np.any(counts > 1):
        ambiguous_query = unique_query_matches[np.where(counts > 1)[0][0]]
        raise InvalidFileError(
            f"The covalent bond in the 'struct_conn' category at index "
            f"{ambiguous_query} cannot be unambiguously assigned to atoms in "
            f"the 'atom_site' category"
        )

    # -1 indicates that no match was found in the reference
    match_indices = np.full(len(query_arrays[0]), -1, dtype=int)
    match_indices[query_matches] = reference_matches
    return match_indices


def _find_matches_by_dict(query_arrays, reference_arrays):
    """
    For each index in the `query_arrays` find the indices in the
    `reference_arrays` where all query values match the reference counterpart.
    If no match is found for a query, the corresponding index is -1.
    """
    # Convert reference arrays to a dictionary for O(1) lookups
    reference_dict = {}
    ambiguous_keys = set()
    for ref_idx, ref_row in enumerate(zip(*reference_arrays)):
        ref_key = tuple(ref_row)
        if ref_key in reference_dict:
            ambiguous_keys.add(ref_key)
            continue
        reference_dict[ref_key] = ref_idx

    match_indices = []
    for query_idx, query_row in enumerate(zip(*query_arrays)):
        query_key = tuple(query_row)
        occurrence = reference_dict.get(query_key)

        if occurrence is None:
            # -1 indicates that no match was found in the reference
            match_indices.append(-1)
        elif query_key in ambiguous_keys:
            # The query cannot be uniquely matched to an atom in the reference
            raise InvalidFileError(
                f"The covalent bond in the 'struct_conn' category at index "
                f"{query_idx} cannot be unambiguously assigned to atoms in "
                f"the 'atom_site' category"
            )
        else:
            match_indices.append(occurrence)

    return np.array(match_indices)


def _find_matches(query_arrays, reference_arrays):
    """
    For each index in the `query_arrays` find the indices in the
    `reference_arrays` where all query values match the reference counterpart.
    If no match is found for a query, the corresponding index is -1.
    """
    # it was observed that when the size exceeds 2^13 (8192)
    # the dict strategy becomes significantly faster than the dense array
    # and does not cause excessive memory usage.
    # https://github.com/biotite-dev/biotite/pull/765#issuecomment-2696591338
    if (
        query_arrays[0].shape[0] * reference_arrays[0].shape[0]
        <= FIND_MATCHES_SWITCH_THRES
    ):
        match_indices = _find_matches_by_dense_array(query_arrays, reference_arrays)
    else:
        match_indices = _find_matches_by_dict(query_arrays, reference_arrays)
    return match_indices


def _parse_inter_residue_bonds(atom_site, struct_conn):
    """
    Create inter-residue bonds by parsing the ``struct_conn`` category.
    The atom indices of each bond are found by matching the bond labels
    to the ``atom_site`` category.
    """
    # Identity symmetry operation
    IDENTITY = "1_555"
    # Columns in 'atom_site' that should be matched by 'struct_conn'
    COLUMNS = [
        "label_asym_id",
        "label_comp_id",
        "label_seq_id",
        "label_atom_id",
        "label_alt_id",
        "auth_asym_id",
        "auth_comp_id",
        "auth_seq_id",
        "pdbx_PDB_ins_code",
    ]

    covale_mask = np.isin(
        struct_conn["conn_type_id"].as_array(str),
        list(PDBX_BOND_TYPE_ID_TO_TYPE.keys()),
    )
    if "ptnr1_symmetry" in struct_conn:
        covale_mask &= struct_conn["ptnr1_symmetry"].as_array(str, IDENTITY) == IDENTITY
    if "ptnr2_symmetry" in struct_conn:
        covale_mask &= struct_conn["ptnr2_symmetry"].as_array(str, IDENTITY) == IDENTITY

    atom_indices = [None] * 2
    for i in range(2):
        reference_arrays = []
        query_arrays = []
        for col_name in COLUMNS:
            struct_conn_col_name = pdbx_convert._get_struct_conn_col_name(
                col_name, i + 1
            )
            if col_name not in atom_site or struct_conn_col_name not in struct_conn:
                continue
            # Ensure both arrays have the same dtype to allow comparison
            reference = atom_site[col_name].as_array()
            dtype = reference.dtype
            query = struct_conn[struct_conn_col_name].as_array(dtype)
            if np.issubdtype(reference.dtype, str):
                # The mask value is not necessarily consistent
                # between query and reference
                # -> make it consistent
                reference[reference == "?"] = "."
                query[query == "?"] = "."
            reference_arrays.append(reference)
            query_arrays.append(query[covale_mask])
        # Match the combination of 'label_asym_id', 'label_comp_id', etc.
        # in 'atom_site' and 'struct_conn'
        atom_indices[i] = _find_matches(query_arrays, reference_arrays)
    atoms_indices_1 = atom_indices[0]
    atoms_indices_2 = atom_indices[1]

    # Some bonds in 'struct_conn' may not be found in 'atom_site'
    # This is okay,
    # as 'atom_site' might already be reduced to a single model
    mapping_exists_mask = (atoms_indices_1 != -1) & (atoms_indices_2 != -1)
    atoms_indices_1 = atoms_indices_1[mapping_exists_mask]
    atoms_indices_2 = atoms_indices_2[mapping_exists_mask]

    bond_type_id = struct_conn["conn_type_id"].as_array()
    # Consecutively apply the same masks as applied to the atom indices
    # Logical combination does not work here,
    # as the second mask was created based on already filtered data
    bond_type_id = bond_type_id[covale_mask][mapping_exists_mask]
    # The type ID is always present in the dictionary,
    # as it was used to filter the applicable bonds
    bond_types = [PDBX_BOND_TYPE_ID_TO_TYPE[type_id] for type_id in bond_type_id]

    return BondList(
        atom_site.row_count,
        np.stack([atoms_indices_1, atoms_indices_2, bond_types], axis=-1),
    )


# copy from biotite 1.1.0
def concatenate(atoms):
    """
    Concatenate multiple :class:`AtomArray` or :class:`AtomArrayStack` objects into
    a single :class:`AtomArray` or :class:`AtomArrayStack`, respectively.

    Parameters
    ----------
    atoms : iterable object of AtomArray or AtomArrayStack
        The atoms to be concatenated.
        :class:`AtomArray` cannot be mixed with :class:`AtomArrayStack`.

    Returns
    -------
    concatenated_atoms : AtomArray or AtomArrayStack
        The concatenated atoms, i.e. its ``array_length()`` is the sum of the
        ``array_length()`` of the input ``atoms``.

    Notes
    -----
    The following rules apply:

    - Only the annotation categories that exist in all elements are transferred.
    - The box of the first element that has a box is transferred, if any.
    - The bonds of all elements are concatenated, if any element has associated bonds.
      For elements without a :class:`BondList` an empty :class:`BondList` is assumed.

    Examples
    --------

    >>> atoms1 = array([
    ...     Atom([1,2,3], res_id=1, atom_name="N"),
    ...     Atom([4,5,6], res_id=1, atom_name="CA"),
    ...     Atom([7,8,9], res_id=1, atom_name="C")
    ... ])
    >>> atoms2 = array([
    ...     Atom([1,2,3], res_id=2, atom_name="N"),
    ...     Atom([4,5,6], res_id=2, atom_name="CA"),
    ...     Atom([7,8,9], res_id=2, atom_name="C")
    ... ])
    >>> print(concatenate([atoms1, atoms2]))
                1      N                1.000    2.000    3.000
                1      CA               4.000    5.000    6.000
                1      C                7.000    8.000    9.000
                2      N                1.000    2.000    3.000
                2      CA               4.000    5.000    6.000
                2      C                7.000    8.000    9.000
    """
    # Ensure that the atoms can be iterated over multiple times
    if not isinstance(atoms, Sequence):
        atoms = list(atoms)

    length = 0
    depth = None
    element_type = None
    common_categories = set(atoms[0].get_annotation_categories())
    box = None
    has_bonds = False
    for element in atoms:
        if element_type is None:
            element_type = type(element)
        else:
            if not isinstance(element, element_type):
                raise TypeError(
                    f"Cannot concatenate '{type(element).__name__}' "
                    f"with '{element_type.__name__}'"
                )
        length += element.array_length()
        if isinstance(element, AtomArrayStack):
            if depth is None:
                depth = element.stack_depth()
            else:
                if element.stack_depth() != depth:
                    raise IndexError("The stack depths are not equal")
        common_categories &= set(element.get_annotation_categories())
        if element.box is not None and box is None:
            box = element.box
        if element.bonds is not None:
            has_bonds = True

    if element_type == AtomArray:
        concat_atoms = AtomArray(length)
    elif element_type == AtomArrayStack:
        concat_atoms = AtomArrayStack(depth, length)
    concat_atoms.coord = np.concatenate([element.coord for element in atoms], axis=-2)
    for category in common_categories:
        concat_atoms.set_annotation(
            category,
            np.concatenate(
                [element.get_annotation(category) for element in atoms], axis=0
            ),
        )
    concat_atoms.box = box
    if has_bonds:
        # Concatenate bonds of all elements
        concat_atoms.bonds = bond_list_concatenate(
            [
                (
                    element.bonds
                    if element.bonds is not None
                    else BondList(element.array_length())
                )
                for element in atoms
            ]
        )

    return concat_atoms


# copy from biotite 1.1.0 and modified
def bond_list_concatenate(bonds_lists):
    """
    Concatenate multiple :class:`BondList` objects into a single
    :class:`BondList`, respectively.

    Parameters
    ----------
    bonds_lists : iterable object of BondList
        The bond lists to be concatenated.

    Returns
    -------
    concatenated_bonds : BondList
        The concatenated bond lists.

    Examples
    --------

    >>> bonds1 = BondList(2, np.array([(0, 1)]))
    >>> bonds2 = BondList(3, np.array([(0, 1), (0, 2)]))
    >>> merged_bonds = BondList.concatenate([bonds1, bonds2])
    >>> print(merged_bonds.get_atom_count())
    5
    >>> print(merged_bonds.as_array()[:, :2])
    [[0 1]
        [2 3]
        [2 4]]
    """
    # Ensure that the bonds_lists can be iterated over multiple times
    if not isinstance(bonds_lists, Sequence):
        bonds_lists = list(bonds_lists)

    merged_bonds = np.concatenate([bond_list._bonds for bond_list in bonds_lists])
    # Offset the indices of appended bonds list
    # (consistent with addition of AtomArray)
    start = 0
    stop = 0
    cum_atom_count = 0
    for bond_list in bonds_lists:
        stop = start + bond_list._bonds.shape[0]
        merged_bonds[start:stop, :2] += cum_atom_count
        cum_atom_count += bond_list._atom_count
        start = stop

    merged_bond_list = BondList(cum_atom_count)
    # Array is not used in constructor to prevent unnecessary
    # maximum and redundant bond calculation
    merged_bond_list._bonds = merged_bonds
    merged_bond_list._max_bonds_per_atom = max(
        [bond_list._max_bonds_per_atom for bond_list in bonds_lists]
    )
    return merged_bond_list


def move_inside_box(coord, box):
    r"""
    Move all coordinates into the given box, with the box vectors
    originating at *(0,0,0)*.

    Coordinates are outside the box, when they cannot be represented by
    a linear combination of the box vectors with scalar factors
    :math:`0 \le a_i \le 1`.
    In this case the affected coordinates are translated by the box
    vectors, so that they are inside the box.

    Parameters
    ----------
    coord : ndarray, dtype=float, shape=(n,3) or shape=(m,n,3)
        The coordinates for one or multiple models.
    box : ndarray, dtype=float, shape=(3,3) or shape=(m,3,3)
        The box(es) for one or multiple models.
        When `coord` is given for multiple models, :attr:`box` must be
        given for multiple models as well.

    Returns
    -------
    moved_coord : ndarray, dtype=float, shape=(n,3) or shape=(m,n,3)
        The moved coordinates.
        Has the same shape is the input `coord`.
    translation_operation_id : int
        The translational offset applied to the operation_ids.

    Examples
    --------

    >>> box = np.array([[10,0,0], [0,10,0], [0,0,10]], dtype=float)
    >>> inside_coord        = [ 1,  2,  3]
    >>> outside_coord       = [ 1, 22, 54]
    >>> other_outside_coord = [-4,  8,  6]
    >>> coord = np.stack([inside_coord, outside_coord, other_outside_coord])
    >>> print(coord)
    [[ 1  2  3]
     [ 1 22 54]
     [-4  8  6]]
    >>> moved_coord, translation_operation_id = move_inside_box(coord, box)
    >>> print(moved_coord.astype(int))
    [[1 2 3]
     [1 2 4]
     [6 8 6]]
    >>> print(translation_operation_id)
    [0 -25 100]
    """
    fractions = coord_to_fraction(coord, box)
    fractions_rem = fractions % 1
    offset = (fractions_rem - fractions).astype(int)
    assert np.all(np.abs(offset) < 5)
    translation_operation_id = offset[:, 0] * 100 + offset[:, 1] * 10 + offset[:, 2]
    return fraction_to_coord(fractions_rem, box), translation_operation_id
