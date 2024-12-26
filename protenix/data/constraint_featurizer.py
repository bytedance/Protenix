import torch
import torch.nn.functional as F
import numpy as np
import pandas as pd
from biotite.structure import Atom, AtomArray, get_residue_starts
from scipy.spatial.distance import cdist

from protenix.data.tokenizer import AtomArrayTokenizer, Token, TokenArray


class ConstraintFeaturizer(object):
    def __init__(
        self,
        token_array: TokenArray,
        atom_array: AtomArray,
        pad_value: float = 0,
        generator=None,
    ):
        self.token_array = token_array
        self.atom_array = atom_array
        self.pad_value = pad_value
        self.generator = generator

        token_centre_atom_indices = self.token_array.get_annotation("centre_atom_index")
        centre_atoms = self.atom_array[token_centre_atom_indices]
        self.asymid = torch.tensor(centre_atoms.asym_id_int, dtype=torch.long)

    @staticmethod
    def one_hot_encoder(feature, num_classes):
        return F.one_hot(feature, num_classes=num_classes).float()

    def encode(self, feature, feature_type, **kwargs):
        if feature_type == "one_hot":
            return ConstraintFeaturizer.one_hot_encoder(
                feature, num_classes=kwargs.get("num_classes", -1)
            )
        elif feature_type == "continuous":
            return feature
        else:
            raise RuntimeError(f"Invalid feature_type: {feature_type}")

    def generate_spec_constraint(
        self,
    ):
        pass


class ContactFeaturizer(ConstraintFeaturizer):
    def __init__(self, **kargs):
        super().__init__(**kargs)

    def generate_spec_constraint(self, contact_specifics, feature_type):
        """
        parse constraint from user specification
        """

        contact_feature = torch.full(
            (self.asymid.shape[0], self.asymid.shape[0], 2),
            fill_value=self.pad_value,
            dtype=torch.float32,
        )
        for token_list_1, token_list_2, max_distance in contact_specifics:

            token_id_1 = token_list_1[
                torch.randint(
                    high=token_list_1.shape[0], size=(1,), generator=self.generator
                ).item()
            ]
            token_id_2 = token_list_2[
                torch.randint(
                    high=token_list_2.shape[0], size=(1,), generator=self.generator
                ).item()
            ]

            contact_feature[token_id_1, token_id_2, 1] = max_distance
            contact_feature[token_id_2, token_id_1, 1] = max_distance
            contact_feature[token_id_1, token_id_2, 0] = 0
            contact_feature[token_id_2, token_id_1, 0] = 0

        contact_feature = self.encode(
            feature=contact_feature, feature_type=feature_type
        )
        return contact_feature
