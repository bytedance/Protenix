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

import unittest
from contextlib import ExitStack
from unittest.mock import Mock, patch

import torch

from protenix.data.pipeline import dataset as dataset_module


class _FakeFeaturizer:
    def __init__(self, **kwargs):
        pass

    def get_all_input_features(self):
        return {}

    def get_labels(self):
        return {}

    def get_atom_permutation_list(self):
        return []

    @staticmethod
    def get_gt_full_complex_features(**kwargs):
        return {}, kwargs["atom_array"]


class TestESMTrainingFeatures(unittest.TestCase):
    def test_dataset_adds_precomputed_esm_features_after_constraints(self):
        dataset = object.__new__(dataset_module.BaseSingleDataset)
        dataset.constraint = {"enable": True}
        dataset.esm_featurizer = Mock(return_value=torch.ones((2, 2560)))
        dataset.find_pocket = False
        dataset.find_eval_chain_interface = False
        dataset.is_distillation = False
        dataset.ref_pos_augment = False
        dataset.lig_atom_rename = False

        token_array = object()
        atom_array = object()
        constrained_token_array = object()
        constrained_atom_array = object()
        full_atom_array = object()
        bioassembly_dict = {"sequences": {"1": "AA"}, "pdb_id": "test"}

        with ExitStack() as stack:
            stack.enter_context(patch.object(dataset_module, "Featurizer", _FakeFeaturizer))
            stack.enter_context(
                patch.object(
                    dataset_module,
                    "make_dummy_feature",
                    side_effect=lambda **kw: kw["features_dict"],
                )
            )
            stack.enter_context(
                patch.object(
                    dataset_module,
                    "data_type_transform",
                    side_effect=lambda **kw: kw["feat_or_label_dict"],
                )
            )
            constraint_feature = stack.enter_context(
                patch.object(
                    dataset,
                    "get_constraint_feature",
                    return_value=(
                        constrained_token_array,
                        constrained_atom_array,
                        {"constraint_feature": {}},
                        {},
                        full_atom_array,
                    ),
                )
            )
            features, _, _ = dataset.get_feature_and_label(
                idx=0,
                token_array=token_array,
                atom_array=atom_array,
                msa_features={},
                template_features={},
                full_atom_array=full_atom_array,
                bioassembly_dict=bioassembly_dict,
            )

        constraint_feature.assert_called_once()
        self.assertIs(features["esm_token_embedding"], dataset.esm_featurizer.return_value)
        dataset.esm_featurizer.assert_called_once_with(
            token_array=constrained_token_array,
            atom_array=constrained_atom_array,
            bioassembly_dict=bioassembly_dict,
            inference_mode=False,
        )

    def test_disabled_esm_does_not_create_a_featurizer(self):
        with patch.object(dataset_module, "ESMFeaturizer") as esm_featurizer:
            result = dataset_module.get_esm_featurizer({"esm": {"enable": False}}, None)
        self.assertIsNone(result)
        esm_featurizer.assert_not_called()

    def test_enabled_esm_uses_the_training_cache_layout(self):
        with patch.object(dataset_module, "ESMFeaturizer") as esm_featurizer:
            dataset_module.get_esm_featurizer(
                {"esm": {"enable": True, "model_name": "esm2-3b"}}, "errors"
            )
        esm_featurizer.assert_called_once_with(
            embedding_dir="./esm_embeddings/esm2-3b",
            sequence_fpath="./scripts/msa/data/pdb_seqs/pdb_labels_seqs.csv",
            embedding_dim=2560,
            error_dir="errors",
        )


if __name__ == "__main__":
    unittest.main()
