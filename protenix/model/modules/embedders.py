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

from typing import Any, Optional, Union

import torch
import torch.nn as nn
import torch.nn.functional as F

from protenix.model.modules.primitives import LinearNoBias
from protenix.model.modules.transformer import AtomAttentionEncoder


class InputFeatureEmbedder(nn.Module):
    """
    Implements Algorithm 2 in AF3
    """

    def __init__(
        self,
        c_atom: int = 128,
        c_atompair: int = 16,
        c_token: int = 384,
        esm_configs: dict = {},
    ) -> None:
        """
        Args:
            c_atom (int, optional): atom embedding dim. Defaults to 128.
            c_atompair (int, optional): atom pair embedding dim. Defaults to 16.
            c_token (int, optional): token embedding dim. Defaults to 384.
        """
        super(InputFeatureEmbedder, self).__init__()
        self.c_atom = c_atom
        self.c_atompair = c_atompair
        self.c_token = c_token
        self.atom_attention_encoder = AtomAttentionEncoder(
            c_atom=c_atom,
            c_atompair=c_atompair,
            c_token=c_token,
            has_coords=False,
        )

        self.esm_configs = {
            "enable": esm_configs.get("enable", False),
            "embedding_dim": esm_configs.get("embedding_dim", 2560),
        }
        if self.esm_configs["enable"]:
            self.linear_esm = LinearNoBias(
                self.esm_configs["embedding_dim"],
                self.c_token + 32 + 32 + 1,
            )
            nn.init.zeros_(self.linear_esm.weight)

        # Line2
        self.input_feature = {"restype": 32, "profile": 32, "deletion_mean": 1}

    def forward(
        self,
        input_feature_dict: dict[str, Any],
        inplace_safe: bool = False,
        chunk_size: Optional[int] = None,
    ) -> torch.Tensor:
        """
        Args:
            input_feature_dict (Dict[str, Any]): dict of input features
            inplace_safe (bool): Whether it is safe to use inplace operations. Defaults to False.
            chunk_size (Optional[int]): Chunk size for memory-efficient operations. Defaults to None.

        Returns:
            torch.Tensor: token embedding
                [..., N_token, 384 (c_token) + 32 + 32 + 1 :=449]
        """
        # Embed per-atom features.
        a, _, _, _ = self.atom_attention_encoder(
            input_feature_dict["atom_to_token_idx"],
            input_feature_dict["ref_pos"],
            input_feature_dict["ref_charge"],
            input_feature_dict["ref_mask"],
            input_feature_dict["ref_atom_name_chars"],
            input_feature_dict["ref_element"],
            input_feature_dict["d_lm"],
            input_feature_dict["v_lm"],
            input_feature_dict["pad_info"],
            inplace_safe=inplace_safe,
            chunk_size=chunk_size,
        )  # [..., N_token, c_token]
        # Concatenate the per-token features.
        batch_shape = input_feature_dict["restype"].shape[:-1]
        s_inputs = torch.cat(
            [a]
            + [
                input_feature_dict[name].reshape(*batch_shape, d)
                for name, d in self.input_feature.items()
            ],
            dim=-1,
        )

        if self.esm_configs["enable"]:
            # Add esm embedding to s_inputs if enable.
            esm_embeddings = self.linear_esm(input_feature_dict["esm_token_embedding"])
            s_inputs = s_inputs + esm_embeddings

        return s_inputs


class RelativePositionEncoding(nn.Module):
    """
    Implements Algorithm 3 in AF3
    """

    def __init__(self, r_max: int = 32, s_max: int = 2, c_z: int = 128) -> None:
        """
        Args:
            r_max (int, optional): Relative position indices clip value. Defaults to 32.
            s_max (int, optional): Relative chain indices clip value. Defaults to 2.
            c_z (int, optional): hidden dim [for pair embedding]. Defaults to 128.
        """
        super(RelativePositionEncoding, self).__init__()
        self.r_max = r_max
        self.s_max = s_max
        self.c_z = c_z
        self.linear_no_bias = LinearNoBias(
            in_features=(4 * self.r_max + 2 * self.s_max + 7), out_features=self.c_z
        )
        self.input_feature = {
            "asym_id": 1,
            "residue_index": 1,
            "entity_id": 1,
            "sym_id": 1,
            "token_index": 1,
        }

    def generate_relp(self, input_feature_dict: dict[str, Any]) -> dict[str, Any]:
        """asym_id, residue_index, entity_id, token_index, sym_id"""
        with torch.no_grad():
            asym_id = input_feature_dict["asym_id"]
            residue_index = input_feature_dict["residue_index"]
            entity_id = input_feature_dict["entity_id"]
            token_index = input_feature_dict["token_index"]
            sym_id = input_feature_dict["sym_id"]

            b_same_chain = (
                asym_id[..., :, None] == asym_id[..., None, :]
            ).long()  # [..., N_token, N_token]
            b_same_residue = (
                residue_index[..., :, None] == residue_index[..., None, :]
            ).long()  # [..., N_token, N_token]
            b_same_entity = (
                entity_id[..., :, None] == entity_id[..., None, :]
            ).long()  # [..., N_token, N_token]
            d_residue = torch.clip(
                input=residue_index[..., :, None]
                - residue_index[..., None, :]
                + self.r_max,
                min=0,
                max=2 * self.r_max,
            ) * b_same_chain + (1 - b_same_chain) * (
                2 * self.r_max + 1
            )  # [..., N_token, N_token]
            a_rel_pos = F.one_hot(d_residue, 2 * (self.r_max + 1))
            d_token = torch.clip(
                input=token_index[..., :, None]
                - token_index[..., None, :]
                + self.r_max,
                min=0,
                max=2 * self.r_max,
            ) * b_same_chain * b_same_residue + (1 - b_same_chain * b_same_residue) * (
                2 * self.r_max + 1
            )  # [..., N_token, N_token]
            a_rel_token = F.one_hot(d_token, 2 * (self.r_max + 1))
            d_chain = torch.clip(
                input=sym_id[..., :, None] - sym_id[..., None, :] + self.s_max,
                min=0,
                max=2 * self.s_max,
            ) * b_same_entity + (1 - b_same_entity) * (
                2 * self.s_max + 1
            )  # [..., N_token, N_token]
            a_rel_chain = F.one_hot(d_chain, 2 * (self.s_max + 1))

            relp = torch.cat(
                [a_rel_pos, a_rel_token, b_same_entity[..., None], a_rel_chain],
                dim=-1,
            ).float()
            input_feature_dict["relp"] = relp
        return input_feature_dict

    def forward(self, relp_feature) -> torch.Tensor:
        """
        Args:
            asym_id / residue_index / entity_id / sym_id / token_index
                [..., N_tokens]
        Returns:
            torch.Tensor: relative position encoding
                [..., N_token, N_token, c_z]
        """
        return self.linear_no_bias(relp_feature)


class FourierEmbedding(nn.Module):
    """
    Implements Algorithm 22 in AF3
    """

    def __init__(self, c: int, seed: int = 42) -> None:
        """
        Args:
            c (int): embedding dim.
        """
        super(FourierEmbedding, self).__init__()
        self.c = c
        self.seed = seed
        generator = torch.Generator()
        generator.manual_seed(seed)
        w_value = torch.randn(size=(c,), generator=generator)
        self.w = nn.Parameter(w_value, requires_grad=False)
        b_value = torch.randn(size=(c,), generator=generator)
        self.b = nn.Parameter(b_value, requires_grad=False)

    def forward(self, t_hat_noise_level: torch.Tensor) -> torch.Tensor:
        """
        Args:
            t_hat_noise_level (torch.Tensor): the noise level
                [..., N_sample]

        Returns:
            torch.Tensor: the output fourier embedding
                [..., N_sample, c]
        """
        return torch.cos(
            input=2 * torch.pi * (t_hat_noise_level.unsqueeze(dim=-1) * self.w + self.b)
        )


class SubstructureEmbedder(nn.Module):
    """
    Implements Substructure Embedder
    """

    def __init__(
        self,
        n_classes: int,
        c_pair_dim: int,
        architecture: str = "mlp",
        hidden_dim: int = 256,
        n_layers: int = 3,
        dropout: float = 0.1,
    ) -> None:
        """
        Args:
            n_classes (int): Number of distance classes in input
            c_pair_dim (int): Output pair embedding dimension
            architecture (str): Either 'mlp' or 'transformer'
            hidden_dim (int): Hidden dimension for both architectures
            n_layers (int): Number of layers (MLP or Transformer)
            dropout (float): Dropout rate
        """
        super().__init__()
        self.architecture = architecture.lower()
        if self.architecture == "mlp":
            layers = []
            layers.append(LinearNoBias(n_classes, hidden_dim))
            layers.append(nn.ReLU())
            layers.append(nn.Dropout(dropout))
            for _ in range(n_layers - 2):
                layers.append(LinearNoBias(hidden_dim, hidden_dim))
                layers.append(nn.ReLU())
                layers.append(nn.Dropout(dropout))

            # Output layer
            layers.append(LinearNoBias(hidden_dim, c_pair_dim))

            self.network = nn.Sequential(*layers)

        elif self.architecture == "transformer":
            self.input_proj = LinearNoBias(n_classes, hidden_dim)
            encoder_layer = nn.TransformerEncoderLayer(
                d_model=hidden_dim,
                nhead=4,
                dim_feedforward=hidden_dim * 4,
                dropout=dropout,
                batch_first=True,
            )
            self.transformer = nn.TransformerEncoder(encoder_layer, num_layers=n_layers)
            self.output_proj = LinearNoBias(hidden_dim, c_pair_dim)
        else:
            raise ValueError(f"Unknown architecture: {architecture}")

    def forward(self, x: torch.Tensor) -> torch.Tensor:
        """
        Args:
            x (torch.Tensor): Input distance map
                shape: [..., N_token, N_token, N_classes]
        Returns:
            torch.Tensor: Output pair embeddings
                shape: [..., N_token, N_token, c_pair_dim]
        """
        if self.architecture == "mlp":
            return self.network(x)
        else:  # transformer
            # Reshape for transformer
            orig_shape = x.shape
            x = x.view(-1, orig_shape[-3], orig_shape[-2], orig_shape[-1])
            batch_size, n_token, _, n_classes = x.shape

            # Project and reshape to sequence
            x = self.input_proj(x)
            x = x.reshape(batch_size, n_token * n_token, -1)

            # Apply transformer
            x = self.transformer(x)

            # Project to output dim and reshape back
            x = self.output_proj(x)
            x = x.reshape(*orig_shape[:-1], -1)

            return x


class ConstraintEmbedder(nn.Module):
    """
    Implements Constraint Embedder
    """

    def __init__(
        self,
        pocket_embedder: dict[str:int],
        contact_embedder: dict[str:int],
        contact_atom_embedder: dict[str:int],
        substructure_embedder: dict[str:int],
        c_constraint_z: int,
        initialize_method: str = "zero",
        **kwarg,
    ) -> None:
        """

        Args:
            pocket_embedder (dict[str:int]): pocket embedder config
            contact_embedder (dict[str:int]): contact embedder config
            contact_atom_embedder (dict[str:int]): contact atom embedder config
            substructure_embedder (dict[str:int]): substructure embedder config
            c_constraint_z (int): constraint z dimension
            initialize_method (str): initialize method
        """
        super(ConstraintEmbedder, self).__init__()
        self.pocket_embedder_config = pocket_embedder
        self.contact_embedder_config = contact_embedder
        self.contact_atom_embedder_config = contact_atom_embedder
        self.substructure_embedder_config = substructure_embedder
        # pocket embedder
        if self.pocket_embedder_config.get("enable", False):
            self.pocket_z_embedder = LinearNoBias(
                in_features=self.pocket_embedder_config.get("c_z_input", 1),
                out_features=c_constraint_z,
            )

        # token contact embedder
        if self.contact_embedder_config.get("enable", False):
            self.contact_z_embedder = LinearNoBias(
                in_features=contact_embedder["c_z_input"], out_features=c_constraint_z
            )

        # atom contact embedder
        if self.contact_atom_embedder_config.get("enable", False):
            self.contact_atom_z_embedder = LinearNoBias(
                in_features=contact_atom_embedder["c_z_input"],
                out_features=c_constraint_z,
            )

        # substructure embedder
        if self.substructure_embedder_config.get("enable", False):
            self.substructure_z_embedder = SubstructureEmbedder(
                n_classes=self.substructure_embedder_config.get("n_classes", 4),
                c_pair_dim=c_constraint_z,
                architecture=self.substructure_embedder_config.get(
                    "architecture", "mlp"
                ),
                hidden_dim=self.substructure_embedder_config.get("hidden_dim", 256),
                n_layers=self.substructure_embedder_config.get("n_layers", 3),
            )

        for module in self.modules():
            if isinstance(module, nn.Linear):
                if initialize_method == "zero":
                    nn.init.zeros_(module.weight)

    def forward(
        self,
        constraint_feature_dict: dict[str, Union[torch.Tensor, int, float, dict]],
    ) -> torch.Tensor:
        """

        Args:
            constraint_feature_dict (dict[str, Union[torch.Tensor, int, float, dict]]): dict of input features

        Returns:
            torch.Tensor: token embedding
                [..., N_token, c_s]
        """
        z_constraint = None

        if self.pocket_embedder_config.get("enable", False):
            z_constraint = self.pocket_z_embedder(constraint_feature_dict["pocket"])

        if self.contact_embedder_config.get("enable", False):
            z_contact = self.contact_z_embedder(constraint_feature_dict["contact"])
            z_constraint = (
                z_contact if z_constraint is None else z_constraint + z_contact
            )

        if self.contact_atom_embedder_config.get("enable", False):
            z_contact_atom = self.contact_atom_z_embedder(
                constraint_feature_dict["contact_atom"]
            )
            z_constraint = (
                z_contact_atom
                if z_constraint is None
                else z_constraint + z_contact_atom
            )

        # substructure embedder
        if self.substructure_embedder_config.get("enable", False):
            z_substructure = self.substructure_z_embedder(
                constraint_feature_dict["substructure"]
            )
            z_constraint = (
                z_substructure
                if z_constraint is None
                else z_constraint + z_substructure
            )
        return z_constraint
