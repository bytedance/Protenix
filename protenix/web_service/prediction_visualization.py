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

import glob
import json
import os

import biotite
import biotite.structure as struc
import matplotlib.pyplot as plt
import numpy as np
import py3Dmol
from biotite.structure.io.pdbx import CIFFile, get_structure
from matplotlib.colors import ListedColormap
from matplotlib.ticker import MaxNLocator


class PredictionLoader:

    def __init__(self, pred_fpath: str):
        self.pred_fpath = pred_fpath
        self._load_cif()
        self._load_confidence_pred()

    def _load_json(self, fpath: str):

        try:
            with open(fpath, "r") as file:
                data = json.load(file)
                return data
        except FileNotFoundError:
            print(f"Error: File '{fpath}' not found.")
        except json.JSONDecodeError:
            print(f"Error: Failed to decode JSON from file '{fpath}'.")
        except Exception as e:
            print(f"An unexpected error occurred: {e}")

    def _convert_to_numpy(self, data: dict):
        for k, v in data.items():
            if isinstance(v, list):
                # convert values to numpy array
                data[k] = np.array(v)
        return data

    def _load_cif(self):
        assert os.path.exists(self.pred_fpath), "prediction file path does not exist."
        fpath_all_preds = glob.glob(os.path.join(self.pred_fpath, "*.cif"))
        self.cif_paths = fpath_all_preds
        self.fnames = [s.split("/")[-1].replace(".cif", "") for s in fpath_all_preds]
        self.preds = [
            get_structure(pdbx_file=CIFFile.read(fpath), model=1, altloc="all")
            for fpath in fpath_all_preds
        ]

    def _load_confidence_pred(self):
        fpath_all_confidences = glob.glob(os.path.join(self.pred_fpath, "*.json"))
        fpath_full_confidences = [
            fpath for fpath in fpath_all_confidences if "full_data_sample" in fpath
        ]
        fpath_summary_confidences = [
            fpath for fpath in fpath_all_confidences if "summary_confidence" in fpath
        ]

        self.summary_confidence_data = [
            self._convert_to_numpy(self._load_json(fpath))
            for fpath in fpath_summary_confidences
        ]

        if fpath_full_confidences:
            self.full_confidence_data = [
                self._convert_to_numpy(self._load_json(fpath))
                for fpath in fpath_full_confidences
            ]
        else:
            self.full_confidence_data = None


def plot_contact_maps_from_pred(
    preds: list,
    fnames: list,
    threshold: float = 7,
    rep_atom: str = "CA",
) -> None:
    """
    Plot contact maps from a directory

    Args
        preds (list): list of biotite.structure.AtomArrays
        fnames (list): list of prediction names
        threshold (int): threshold cutoff for the contact map
        rep_atom (str): type of representative atoms, commonly options are 'CA' and 'CB'
    """

    adjacency_matrices = []
    for pred in preds:
        rep_atom_coord = pred[pred.atom_name == rep_atom]
        cell_list = struc.CellList(rep_atom_coord, cell_size=threshold)
        adjacency_matrices.append(cell_list.create_adjacency_matrix(threshold))

    cmap = ListedColormap(["white", biotite.colors["dimgreen"]])
    fig, axes = plt.subplots(nrows=1, ncols=len(fnames), figsize=(len(fnames) * 3, 6))

    for i in range(len(fnames)):
        axes[i].matshow(adjacency_matrices[i], cmap=cmap, origin="lower")
        axes[i].xaxis.tick_bottom()
        axes[i].set_aspect("equal")
        axes[i].set_xlabel("Residue Number")
        axes[i].set_ylabel("Residue Number")
        axes[i].set_title(fnames[i])

    fig.tight_layout()
    plt.show()


def plot_confidence_measures_from_pred(
    full_confidence: list,
    summary_confidence: list,
    fnames: list,
    show_global_confidence: bool,
    *args,
    **kwargs,
) -> None:
    """
    Plot contact maps from a directory

    Args
        full_confidence (list): list of full confidence metrics from prediction
        summary_confidence (list): list of summary confidence metrics from prediction
        fnames (list): list of prediction names
    """

    atom_plddts = [d["atom_plddt"] for d in full_confidence]
    token_pdes = [d["token_pair_pde"] for d in full_confidence]
    token_paes = [d["token_pair_pae"] for d in full_confidence]
    summary_keys = ["plddt", "gpde", "ptm", "iptm"]

    fig, axes = plt.subplots(
        nrows=3,
        ncols=len(fnames),
        figsize=(len(fnames) * 5, 10),
        gridspec_kw={"height_ratios": [1.5, 2, 2]},
    )
    for i in range(3):
        for j in range(len(fnames)):
            summary_text = ", ".join(
                [
                    f"{k}:{v:.4f}"
                    for k, v in summary_confidence[j].items()
                    if k in summary_keys
                ]
            )
            if i == 0:
                axes[i, j].plot(atom_plddts[j], color="k")
                axes[i, j].set_title(fnames[j], fontsize=15, pad=20)
                if show_global_confidence:
                    axes[i, j].text(
                        0.5,
                        1.07,
                        summary_text,
                        ha="center",
                        va="center",
                        transform=axes[i, j].transAxes,
                        fontsize=10,
                    )
                axes[i, j].set_xlabel("Atom ID", fontsize=12)
                axes[i, j].set_ylabel("pLDDT", fontsize=12)
                axes[i, j].set_ylim([0, 100])
                axes[i, j].spines[["right", "top"]].set_visible(False)

            else:
                data_to_plot = token_pdes[j] if i == 1 else token_paes[j]
                cax = axes[i, j].matshow(data_to_plot, origin="lower")
                axes[i, j].xaxis.tick_bottom()
                axes[i, j].set_aspect("equal")
                axes[i, j].set_xlabel("Scored Residue", fontsize=12)
                axes[i, j].set_ylabel("Aligned Residue", fontsize=12)
                axes[i, j].xaxis.set_major_locator(
                    MaxNLocator(3)
                )  # Max 5 ticks on the x-axis
                axes[i, j].yaxis.set_major_locator(
                    MaxNLocator(3)
                )  # Max 5 ticks on the y-axis
                color_bar = fig.colorbar(
                    cax, ax=axes[i, j], orientation="vertical", pad=0.1, shrink=0.6
                )
                cbar_label = (
                    "Predicted Distance Error" if i == 1 else "Predicted Aligned Error"
                )
                color_bar.set_label(cbar_label)

    fig.tight_layout()
    fig.subplots_adjust(hspace=0.4)

    plt.show()


def plot_3d(
    pred_id=0,
    pred_loader=None,
    show_sidechains=False,
    show_mainchains=False,
    color="rainbow",
):
    view = py3Dmol.view(
        js="https://3dmol.org/build/3Dmol.js",
    )
    fpath = pred_loader.cif_paths[pred_id]
    view.addModelsAsFrames(open(fpath, "r").read(), "cif")
    if color == "pLDDT":
        assert pred_loader is not None
        plddt = pred_loader.full_confidence_data[0]["atom_plddt"]
        for i, score in enumerate(plddt):
            normalized_color = int(
                255 * (score - min(plddt)) / (max(plddt) - min(plddt))
            )
            color = f"rgb({normalized_color}, 0, {255 - normalized_color})"  # Gradient from blue to red

            # Apply color to each atom individually based on pLDDT score
            view.setStyle(
                {"serial": i + 1}, {"cartoon": {"color": color, "min": 50, "max": 90}}
            )

    elif color == "rainbow":
        view.setStyle({"cartoon": {"color": "spectrum"}})

    if show_sidechains:
        BB = ["C", "O", "N"]
        view.addStyle(
            {
                "and": [
                    {"resn": ["GLY", "PRO"], "invert": True},
                    {"atom": BB, "invert": True},
                ]
            },
            {"stick": {"colorscheme": f"WhiteCarbon", "radius": 0.3}},
        )
        view.addStyle(
            {"and": [{"resn": "GLY"}, {"atom": "CA"}]},
            {"sphere": {"colorscheme": f"WhiteCarbon", "radius": 0.3}},
        )
        view.addStyle(
            {"and": [{"resn": "PRO"}, {"atom": ["C", "O"], "invert": True}]},
            {"stick": {"colorscheme": f"WhiteCarbon", "radius": 0.3}},
        )
    if show_mainchains:
        BB = ["C", "O", "N", "CA"]
        view.addStyle(
            {"atom": BB}, {"stick": {"colorscheme": f"WhiteCarbon", "radius": 0.3}}
        )

    view.zoomTo()
    return view
