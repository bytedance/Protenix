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

import copy
import json

import numpy as np
import pytest
from protenix.data.template.template_utils import TemplateHitFeaturizer


@pytest.fixture(scope="module")
def json_template_input():
    with open("examples/example_with_json_template/demo_ab.json") as f:
        query_sequence = json.load(f)[0]["sequences"][0]["proteinChain"]["sequence"]
    with open("examples/example_with_json_template/h_seq.json") as f:
        template = json.load(f)[0]
    return query_sequence, template


@pytest.fixture
def featurizer():
    return TemplateHitFeaturizer(
        mmcif_dir="/tmp/dummy_mmcif_dir",
        template_cache_dir=None,
        kalign_binary_path=None,
        _zero_center_positions=True,
    )


def _parse_mapping(featurizer, json_template_input, query_idx, template_idx):
    query_sequence, source_template = json_template_input
    template = copy.deepcopy(source_template)
    template["queryIndices"] = [query_idx]
    template["templateIndices"] = [template_idx]
    return featurizer.parse_json_templates([template], query_sequence)


def test_json_template_rejects_negative_query_index(featurizer, json_template_input):
    result = _parse_mapping(featurizer, json_template_input, -1, 0)

    assert result.errors == ["Invalid queryIndices at index 0: negative index"]
    assert result.features == []
    assert result.hits == []


def test_json_template_rejects_template_index_below_gap(
    featurizer, json_template_input
):
    result = _parse_mapping(featurizer, json_template_input, 0, -2)

    assert result.errors == ["Invalid templateIndices at index 0: index below -1"]
    assert result.features == []
    assert result.hits == []


def test_json_template_accepts_gap_index(featurizer, json_template_input):
    result = _parse_mapping(featurizer, json_template_input, 0, -1)

    assert result.errors == []
    assert np.sum(result.features[0]["template_all_atom_masks"][0]) == 0
    assert result.features[0]["template_sequence"][:1] == b"-"
    assert result.hits[0].indices_hit[0] == -1


def test_json_template_preserves_normal_mapping(featurizer, json_template_input):
    result = _parse_mapping(featurizer, json_template_input, 0, 0)

    assert result.errors == []
    assert np.sum(result.features[0]["template_all_atom_masks"][0]) > 0
    assert result.features[0]["template_sequence"][:1] == b"Q"
    assert result.hits[0].indices_hit[0] == 0
