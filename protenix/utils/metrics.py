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
import torch

from protenix.utils.distributed import gather_and_merge

common_aggregator = {
    "avg": lambda x: np.mean(x),
    "median": lambda x: np.median(x),
    "pct90": lambda x: np.percentile(x, 90),
    "pct99": lambda x: np.percentile(x, 99),
    "max": lambda x: np.max(x),
    "min": lambda x: np.min(x),
}


class SimpleMetricAggregator(object):
    """A quite simple metrics calculator that only do simple metrics aggregation."""

    def __init__(
        self, aggregator_names=None, gather_before_calc=True, need_gather=True
    ):
        super(SimpleMetricAggregator, self).__init__()
        self.gather_before_calc = gather_before_calc
        self.need_gather = need_gather
        self._metric_data = {}

        self.aggregators = {name: common_aggregator[name] for name in aggregator_names}

    def add(self, key, value, namespace="default"):
        value_dict = self._metric_data.setdefault(namespace, {})
        value_dict.setdefault(key, [])
        if isinstance(value, (float, int)):
            value = np.array([value])
        elif isinstance(value, torch.Tensor):
            if value.dim() == 0:
                value = np.array([value.item()])
            else:
                value = value.detach().cpu().numpy()
        elif isinstance(value, np.ndarray):
            pass
        else:
            raise ValueError(f"Unsupported type for metric data: {type(value)}")
        value_dict[key].append(value)

    def calc(self):
        metric_data, self._metric_data = self._metric_data, {}
        if self.need_gather and self.gather_before_calc:
            metric_data = gather_and_merge(
                metric_data, aggregation_func=lambda l: sum(l, [])
            )
        results = {}
        for agg_name, agg_func in self.aggregators.items():
            for namespace, value_dict in metric_data.items():
                for key, data in value_dict.items():
                    plain_key = f"{namespace}/{key}" if namespace != "default" else key
                    plain_key = f"{plain_key}.{agg_name}"
                    results[plain_key] = agg_func(np.concatenate(data, axis=0))
        if self.need_gather and not self.gather_before_calc:  # need gather after calc
            results = gather_and_merge(results, aggregation_func=np.mean)
        return results
