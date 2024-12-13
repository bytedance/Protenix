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

import re

import ipywidgets as widgets


class DnaRnaProteinEntityWidget(widgets.VBox):
    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.layout.border = "1px solid black"
        self.layout.padding = "1px"
        self.layout.margin = "1px"
        self.molecule_type_dropdown = widgets.Dropdown(
            options=["Protein", "DNA", "RNA"],
            value="Protein",
            disabled=False,
        )
        self.recode_name = {
            "Protein": "proteinChain",
            "DNA": "dnaSequence",
            "RNA": "rnaSequence",
        }

        self.copies_int_text = widgets.BoundedIntText(
            value=1, min=1, max=100, step=1, description="Copies:", disabled=False
        )

        self.sequence_text = widgets.Textarea(
            value="",
            placeholder="Paste sequence or fasta",
            description="Sequence:",
            disabled=False,
        )

        def handle_submit(sender):
            item = sender["owner"]
            cleaned_value = re.sub(r"[^a-zA-Z]", "", item.value)
            item.value = cleaned_value.upper()

        self.sequence_text.observe(handle_submit, names="value")

        self.add_modification = widgets.Button(
            description="modification",
            disabled=False,
            tooltip="add modification",
            icon="plus",
        )

        self.modifications = list()
        self.modifications_id = list()
        self.modifications_count = 0

        self.add_modification.on_click(self.add_modification_callback)

        self.update()

    def add_modification_callback(self, b):
        self.modifications_count += 1
        index = self.modifications_count
        modificationType = widgets.Text(
            value="CCD_XXX",
            placeholder="CCD_XXX",
            description="modificationType",
            max_length=6,
            disabled=False,
        )

        def handle_submit(sender):
            item = sender["owner"]
            item.value = item.value.upper()
            if item.value[0:4] != "CCD_":
                item.value = "CCD_"
            item.value = item.value[0:10]
            if len(item.value) < 4:
                item.value = "CCD_"
            item.value = re.sub(r"[^A-Z0-9_]", "", item.value)

        modificationType.observe(handle_submit, names="value")

        ptmPosition = widgets.BoundedIntText(
            value=1, min=1, description="Position:", disabled=False
        )
        delete_button = widgets.Button(
            description="modification",
            disabled=False,
            button_style="danger",
            tooltip="Click to delete this row",
            icon="trash",
        )
        delete_button.on_click(lambda b: self.delete_ptm_row(index))
        self.modifications.append(
            widgets.HBox([modificationType, ptmPosition, delete_button])
        )
        self.modifications_id.append(index)
        self.update()

    def update(self):
        self.children = (
            [
                widgets.HBox(
                    [
                        self.molecule_type_dropdown,
                        self.copies_int_text,
                        self.sequence_text,
                    ]
                )
            ]
            + [i for i in self.modifications]
            + [self.add_modification]
        )

    def delete_ptm_row(self, index):
        for idx in range(len(self.modifications_id)):
            if self.modifications_id[idx] == index:
                self.modifications.pop(idx)
                self.modifications_id.pop(idx)
                break
        self.update()

    def get_result(self):
        result = dict()
        sequence_type = self.molecule_type_dropdown.value
        recode_name = self.recode_name[sequence_type]
        sequence = self.sequence_text.value
        result[recode_name] = dict()
        result[recode_name]["count"] = self.copies_int_text.value
        result[recode_name]["sequence"] = sequence
        assert len(result[recode_name]["sequence"]) > 0, "sequence length must > 0"
        result[recode_name]["modifications"] = [
            {
                "modificationType": item.children[0].value,
                "Position": item.children[1].value,
            }
            for item in self.modifications
        ]
        for item in result[recode_name]["modifications"]:
            pos = item["Position"]
            assert pos <= len(result[recode_name]["sequence"]), "position out of range"

        if sequence_type == "Protein":
            valid_chars = set("ARNDCQEGHILKMFPSTWYVX")
            for char in sequence:
                assert (
                    char in valid_chars
                ), f"Invalid character '{char}' in sequence. valid set is {valid_chars}"
        elif sequence_type == "DNA":
            valid_chars = set("AGCTXINU")
            for char in sequence:
                assert (
                    char in valid_chars
                ), f"Invalid character '{char}' in sequence. valid set is {valid_chars}"
        elif sequence_type == "RNA":
            valid_chars = set("AGCUXIN")
            for char in sequence:
                assert (
                    char in valid_chars
                ), f"Invalid character '{char}' in sequence. valid set is {valid_chars}"
        return result


class LigandSmilesEntityWidget(widgets.VBox):
    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.layout.border = "1px solid black"
        self.layout.padding = "1px"
        self.layout.margin = "1px"

        self.molecule_type_dropdown = widgets.Dropdown(
            options=["Ligand"],
            value="Ligand",
            disabled=False,
        )

        self.copies_int_text = widgets.BoundedIntText(
            value=1, min=1, max=100, step=1, description="Copies:", disabled=False
        )

        self.smiles_text = widgets.Textarea(
            value="",
            placeholder="SMILES",
            description="SMILES:",
            disabled=False,
        )

        self.update()

    def update(self):
        self.children = [
            widgets.HBox(
                [self.molecule_type_dropdown, self.copies_int_text, self.smiles_text]
            )
        ]

    def get_result(self):
        result = dict()
        name = "ligand"
        result[name] = dict()
        result[name][name] = self.smiles_text.value
        assert len(result[name][name]) > 0, "SMILES length must > 0"
        result[name]["count"] = self.copies_int_text.value
        return result


class LigandIonCCDEntityWidget(widgets.VBox):
    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.layout.border = "1px solid black"
        self.layout.padding = "1px"
        self.layout.margin = "1px"

        self.molecule_type_dropdown = widgets.Dropdown(
            options=["Ligand", "Ion"],
            value="Ligand",
            disabled=False,
        )

        self.copies_int_text = widgets.BoundedIntText(
            value=1, min=1, max=100, step=1, description="Copies:", disabled=False
        )

        self.sequence_text = widgets.Text(
            value="CCD_XXX",
            placeholder="Paste sequence or fasta",
            description="CCD:",
            disabled=False,
        )

        def handle_submit(sender):
            item = sender["owner"]
            item.value = item.value.upper()
            if item.value[0:4] != "CCD_":
                item.value = "CCD_"
            if len(item.value) < 4:
                item.value = "CCD_"

        self.sequence_text.observe(handle_submit, names="value")

        self.update()

    def update(self):
        self.children = [
            widgets.HBox(
                [self.molecule_type_dropdown, self.copies_int_text, self.sequence_text]
            )
        ]

    def get_result(self):
        result = dict()
        if self.molecule_type_dropdown.value == "Ion":
            name = "ion"
        else:
            name = "ligand"
        result[name] = dict()
        result[name][name] = (
            self.sequence_text.value
            if name == "ligand"
            else self.sequence_text.value[4::]
        )
        result[name]["count"] = self.copies_int_text.value
        return result


class CovalentBondsWidget(widgets.VBox):
    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.layout.border = "1px solid black"
        self.layout.padding = "1px"
        self.layout.margin = "1px"

        self.convalent_bonds = list()
        self.convalent_bonds_id = list()
        self.convalent_bonds_count = 0

        self.add_convalent_bond_buttom = widgets.Button(
            description="convalent_bond",
            disabled=False,
            tooltip="add convalent_bond",
            icon="plus",
        )

        self.add_convalent_bond_buttom.on_click(self.on_add_convalent_bond)

        self.update()

    def on_add_convalent_bond(self, b):
        self.convalent_bonds_count += 1
        index = self.convalent_bonds_count

        left_entity = widgets.BoundedIntText(
            value=1, min=1, description="left entity", disabled=False
        )
        left_position = widgets.BoundedIntText(
            value=1, min=1, description="left pos", disabled=False
        )
        left_atom = widgets.Text(
            value="", placeholder="left atom", description="left atom", disabled=False
        )
        right_entity = widgets.BoundedIntText(
            value=1, min=1, description="right entity", disabled=False
        )
        right_position = widgets.BoundedIntText(
            value=1, min=1, description="right pos", disabled=False
        )
        right_atom = widgets.Text(
            value="", placeholder="right atom", description="right atom", disabled=False
        )
        delete_button = widgets.Button(
            description="convalent bond",
            disabled=False,
            button_style="danger",
            tooltip="Click to delete this row",
            icon="trash",
        )
        delete_button.on_click(lambda b: self.delete_row(index))
        self.convalent_bonds.append(
            widgets.HBox(
                [
                    left_entity,
                    left_position,
                    left_atom,
                    right_entity,
                    right_position,
                    right_atom,
                    delete_button,
                ]
            )
        )
        self.convalent_bonds_id.append(index)
        self.update()

    def delete_row(self, index):
        for idx in range(len(self.convalent_bonds_id)):
            if self.convalent_bonds_id[idx] == index:
                self.convalent_bonds.pop(idx)
                self.convalent_bonds_id.pop(idx)
                break
        self.update()

    def update(self):
        self.children = [item for item in self.convalent_bonds] + [
            self.add_convalent_bond_buttom
        ]

    def get_result(self):
        result = [
            {
                "left_entity": item.children[0].value,
                "left_position": item.children[1].value,
                "left_atom": item.children[2].value,
                "right_entity": item.children[3].value,
                "right_position": item.children[4].value,
                "right_atom": item.children[5].value,
            }
            for item in self.convalent_bonds
        ]
        return result


class ModelWidget(widgets.VBox):
    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.layout.border = "1px solid black"
        self.layout.padding = "1px"
        self.layout.margin = "1px"

        self.sample = widgets.BoundedIntText(
            value=5,
            description="N_sample",
            tooltip="diffusion model samples",
            disabled=False,
        )

        self.step = widgets.BoundedIntText(
            value=200,
            max=500,
            description="N_step",
            tooltip="diffusion model step",
            disabled=False,
        )

        self.cycle = widgets.BoundedIntText(
            value=10, description="N_cycle", tooltip="n cycle of model", disabled=False
        )

        self.seeds = widgets.Text(
            value="1,2,3,4,5",
            description="seeds",
            tooltip="random seeds",
            disabled=False,
        )

        self.version = widgets.Dropdown(
            options=["v1", "v2", "v3", "v4", "v5"],
            value="v1",
            description="version",
            tooltip="model version",
        )

        self.update()

    def update(self):
        self.children = [
            widgets.HBox([self.sample, self.step, self.cycle, self.seeds, self.version])
        ]

    def get_result(self):
        results = dict()
        results["model_seeds"] = [int(i) for i in self.seeds.value.split(",")]
        results["model_version"] = self.version.value
        results["N_sample"] = self.sample.value
        results["N_step"] = self.step.value
        results["N_cycle"] = self.cycle.value
        return results


class ProtenixInputViewer(widgets.VBox):
    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.layout.border = "1px solid black"
        self.layout.padding = "1px"
        self.layout.margin = "1px"

        self.name = widgets.Text(
            value="", placeholder="name", description="name", disabled=False
        )

        def handle_submit(sender):
            item = sender["owner"]
            cleaned_value = re.sub(r"[^a-zA-Z0-9_]", "", item.value)
            item.value = cleaned_value

        self.name.observe(handle_submit, names="value")

        self.use_msa = widgets.Checkbox(
            value=True, description="use_msa", disabled=False
        )
        self.atom_confidence = widgets.Checkbox(
            value=False, description="atom_confidence", disabled=False
        )

        self.ligand_smiles_entities = list()
        self.ligand_smiles_entities_id = list()
        self.ligand_smiles_entities_count = 0

        self.dna_rna_protein_entities = list()
        self.dna_rna_protein_entities_id = list()
        self.dna_rna_protein_entities_count = 0

        self.covalent_bonds = CovalentBondsWidget()

        self.add_ligand_smiles = widgets.Button(
            description="Ligand(SMILES)",
            disabled=False,
            tooltip="add Ligand(SMILES) Entity",
            icon="plus",
        )

        self.add_ligand_ion_ccd = widgets.Button(
            description="Ligand/Ion CCD",
            disabled=False,
            tooltip="add Ligand/Ion Entity",
            icon="plus",
        )

        self.add_dna_rna_protein = widgets.Button(
            description="Dna/Rna/Protein",
            disabled=False,
            tooltip="add Dna/Rna/Protein Entity",
            icon="plus",
        )

        self.model_parameter = ModelWidget()

        self.add_ligand_smiles.on_click(self.add_ligand_smiles_callback)
        self.add_ligand_ion_ccd.on_click(self.add_ligand_ion_ccd_callback)
        self.add_dna_rna_protein.on_click(self.add_dna_rna_protein_callback)

        self.update()

    def add_ligand_smiles_callback(self, b):
        self.ligand_smiles_entities_count += 1
        index = self.ligand_smiles_entities_count
        item = LigandSmilesEntityWidget()
        delete_button = widgets.Button(
            description="entity",
            disabled=False,
            button_style="danger",
            tooltip="Click to delete this row",
            icon="trash",
        )
        delete_button.on_click(lambda b: self.delete_ligand_ion_ccd_row(index))
        self.ligand_smiles_entities.append(widgets.HBox([item, delete_button]))
        self.ligand_smiles_entities_id.append(index)
        self.update()

    def add_ligand_ion_ccd_callback(self, b):
        self.ligand_smiles_entities_count += 1
        index = self.ligand_smiles_entities_count
        item = LigandIonCCDEntityWidget()
        delete_button = widgets.Button(
            description="entity",
            disabled=False,
            button_style="danger",
            tooltip="Click to delete this row",
            icon="trash",
        )
        delete_button.on_click(lambda b: self.delete_ligand_ion_ccd_row(index))
        self.ligand_smiles_entities.append(widgets.HBox([item, delete_button]))
        self.ligand_smiles_entities_id.append(index)
        self.update()

    def delete_ligand_ion_ccd_row(self, index):
        for idx in range(len(self.ligand_smiles_entities_id)):
            if self.ligand_smiles_entities_id[idx] == index:
                self.ligand_smiles_entities.pop(idx)
                self.ligand_smiles_entities_id.pop(idx)
                break
        self.update()

    def add_dna_rna_protein_callback(self, b):
        self.dna_rna_protein_entities_count += 1
        index = self.dna_rna_protein_entities_count
        item = DnaRnaProteinEntityWidget()
        delete_button = widgets.Button(
            description="entity",
            disabled=False,
            button_style="danger",
            tooltip="Click to delete this row",
            icon="trash",
        )
        delete_button.on_click(lambda b: self.delete_dna_rna_protein_row(index))
        self.dna_rna_protein_entities.append(widgets.HBox([item, delete_button]))
        self.dna_rna_protein_entities_id.append(index)
        self.update()

    def delete_dna_rna_protein_row(self, index):
        for idx in range(len(self.dna_rna_protein_entities_id)):
            if self.dna_rna_protein_entities_id[idx] == index:
                self.dna_rna_protein_entities.pop(idx)
                self.dna_rna_protein_entities_id.pop(idx)
                break
        self.update()

    def update(self):
        self.children = (
            [widgets.HBox([self.name, self.use_msa, self.atom_confidence])]
            + [item for item in self.dna_rna_protein_entities]
            + [item for item in self.ligand_smiles_entities]
            + [self.covalent_bonds]
            + [
                widgets.HBox(
                    [
                        self.add_dna_rna_protein,
                        self.add_ligand_smiles,
                        self.add_ligand_ion_ccd,
                    ]
                )
            ]
            + [self.model_parameter]
        )

    def get_result(self):
        result = dict()
        result["name"] = self.name.value
        assert len(result["name"]) != 0, "name is empty"
        result["use_msa"] = self.use_msa.value
        result["atom_confidence"] = self.atom_confidence.value
        result["sequences"] = [
            item.children[0].get_result() for item in self.dna_rna_protein_entities
        ] + [item.children[0].get_result() for item in self.ligand_smiles_entities]
        result["covalent_bonds"] = self.covalent_bonds.get_result()
        model_paramter = self.model_parameter.get_result()
        for key in model_paramter:
            result[key] = model_paramter[key]

        entity_num = len(result["sequences"])
        for covalent_bond in result["covalent_bonds"]:
            assert (
                covalent_bond["left_entity"] <= entity_num
                and covalent_bond["left_entity"] > 0
            ), "covalent bond index out of range"
            assert (
                covalent_bond["right_entity"] <= entity_num
                and covalent_bond["right_entity"] > 0
            ), "covalent bond index out of range"

            assert len(covalent_bond["left_atom"]) != 0, "left atom is empty"
            assert len(covalent_bond["right_atom"]) != 0, "right_atom atom is empty"

            left_entiy = result["sequences"][covalent_bond["left_entity"] - 1]
            left_entiy_seq = left_entiy[[item for item in left_entiy.keys()][0]]
            right_entity = result["sequences"][covalent_bond["right_entity"] - 1]
            right_entity_seq = right_entity[[item for item in right_entity.keys()][0]]
            assert (
                covalent_bond["left_position"] <= len(left_entiy_seq)
                and covalent_bond["left_position"] >= 1
            ), "left position out of range"
            assert (
                covalent_bond["right_position"] <= len(right_entity_seq)
                and covalent_bond["right_position"] >= 1
            ), "right position out of range"
        return result
