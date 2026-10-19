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
import torch
import torch.nn.functional as F

from protenix.utils.metrics import SimpleMetricAggregator


@pytest.mark.parametrize("numpy_dtype", [np.float32, np.float64, np.int64])
@pytest.mark.parametrize("gather_before_calc", [False, True])
@pytest.mark.parametrize("need_gather", [False, True])
def test_scalar_numpy_metrics_can_mix_with_other_supported_values(
    numpy_dtype, gather_before_calc, need_gather
):
    aggregator = SimpleMetricAggregator(
        ["avg", "median", "pct90", "pct99", "min", "max"],
        gather_before_calc=gather_before_calc,
        need_gather=need_gather,
    )
    scalar = np.array(1, dtype=numpy_dtype)
    for value in [
        scalar,
        np.array([3, 8], dtype=numpy_dtype),
        np.array(6, dtype=numpy_dtype),
        10,
        torch.tensor(2.0),
    ]:
        aggregator.add("loss", value, namespace="eval")

    assert aggregator.calc() == pytest.approx(
        {
            "eval/loss.avg": 5,
            "eval/loss.median": 4.5,
            "eval/loss.pct90": 9,
            "eval/loss.pct99": 9.9,
            "eval/loss.min": 1,
            "eval/loss.max": 10,
        }
    )
    assert scalar.shape == ()
    assert scalar.dtype == numpy_dtype
    assert aggregator.calc() == {}


@pytest.mark.parametrize("dtype", [torch.float32, torch.float64])
def test_numpy_converted_loss_scalars_match_tensor_metrics(dtype):
    tensor_metrics = SimpleMetricAggregator(["avg"])
    numpy_metrics = SimpleMetricAggregator(["avg"])
    for prediction in [2.0, -1.0, 3.0]:
        loss = F.mse_loss(
            torch.tensor([prediction], dtype=dtype), torch.zeros(1, dtype=dtype)
        )
        tensor_metrics.add("loss", loss)
        numpy_metrics.add("loss", loss.detach().cpu().numpy())

    expected = {"loss.avg": 14 / 3}
    assert tensor_metrics.calc() == pytest.approx(expected)
    assert numpy_metrics.calc() == pytest.approx(expected)


@pytest.mark.parametrize("shape", [(4,), (2, 2)])
@pytest.mark.parametrize("gather_before_calc", [False, True])
@pytest.mark.parametrize("need_gather", [False, True])
def test_batched_metrics_keep_existing_aggregation_and_reset(
    shape, gather_before_calc, need_gather
):
    aggregator = SimpleMetricAggregator(
        ["avg", "max"],
        gather_before_calc=gather_before_calc,
        need_gather=need_gather,
    )
    values = np.array([1.0, 3.0, 7.0, 9.0]).reshape(shape)
    aggregator.add("numpy", values)
    aggregator.add("tensor", torch.from_numpy(values))

    assert aggregator.calc() == pytest.approx(
        {"numpy.avg": 5, "numpy.max": 9, "tensor.avg": 5, "tensor.max": 9}
    )
    aggregator.add("next", np.array([2.0, 4.0]), namespace="train")
    assert aggregator.calc() == pytest.approx(
        {"train/next.avg": 3, "train/next.max": 4}
    )
    assert aggregator.calc() == {}
