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

"""Conservative reconstruction of unobserved branches at bonded oxygen atoms."""

import networkx as nx
import numpy as np
from biotite.structure import BondType


def absent_oxygen_branch(atom_array, component, central_name, central_indices):
    """Return a branch only when every observed copy supports its removal.

    A deposited external single bond can replace one of the two CCD single
    bonds at a neutral oxygen. Only a wholly absent, acyclic side is eligible;
    the retained side must contain deposited atoms in the same CCD component.
    """
    if component.bonds is None or atom_array.bonds is None or not central_indices:
        return []
    if "charge" not in component.get_annotation_categories():
        return []
    matches = np.flatnonzero(component.atom_name == central_name)
    if len(matches) != 1 or len(set(component.atom_name)) != len(component):
        return []
    center = int(matches[0])
    if component.element[center] != "O" or component.charge[center] != 0:
        return []
    neighbors, orders = component.bonds.get_bonds(center)
    heavy = np.isin(component.element[neighbors], ["H", "D"], invert=True)
    neighbors, orders = neighbors[heavy], orders[heavy]
    if len(neighbors) != 2 or np.any(orders != BondType.SINGLE):
        return []
    graph = nx.Graph()
    graph.add_nodes_from(np.flatnonzero(~np.isin(component.element, ["H", "D"])))
    graph.add_edges_from(
        (int(i), int(j))
        for i, j, _ in component.bonds.as_array()
        if i in graph and j in graph
    )
    connected = nx.node_connected_component(graph, center)
    graph.remove_node(center)
    sides = [set(nx.node_connected_component(graph, int(n))) for n in neighbors]
    if sides[0] & sides[1]:  # The oxygen is in a ring.
        return []
    agreed = None
    for idx in central_indices:
        if atom_array.res_name[idx] != component.res_name[center]:
            return []
        same_residue = (
            (atom_array.chain_id == atom_array.chain_id[idx])
            & (atom_array.res_id == atom_array.res_id[idx])
            & (atom_array.ins_code == atom_array.ins_code[idx])
            & (atom_array.res_name == atom_array.res_name[idx])
        )
        bonded, bond_orders = atom_array.bonds.get_bonds(idx)
        external = ~same_residue[bonded]
        partners, external_orders = bonded[external], bond_orders[external]
        if len(partners) != 1 or external_orders[0] != BondType.SINGLE:
            return []
        # Coordination or an unspecified bond cannot certify covalent valence.
        if atom_array.element[partners[0]] not in {"C", "N", "O", "P", "S", "SE"}:
            return []
        observed = set(atom_array.atom_name[same_residue])
        eligible = []
        for side in sides:
            names = set(component.atom_name[list(side)])
            retained = connected - side - {center}
            if (
                not names & observed
                and nx.is_tree(graph.subgraph(side))
                and set(component.atom_name[list(retained)]) & observed
            ):
                eligible.append(names)
        if len(eligible) != 1:
            return []
        if agreed is not None and agreed != eligible[0]:
            return []
        agreed = eligible[0]
    if agreed and "label_entity_id" in atom_array.get_annotation_categories():
        # An entity template is shared even with copies missing the center itself.
        first = central_indices[0]
        same_site = (
            (atom_array.label_entity_id == atom_array.label_entity_id[first])
            & (atom_array.res_id == atom_array.res_id[first])
            & (atom_array.ins_code == atom_array.ins_code[first])
            & (atom_array.res_name == atom_array.res_name[first])
        )
        if set(atom_array.atom_name[same_site]) & agreed:
            return []
    return sorted(agreed) if agreed else []
