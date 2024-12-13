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

import os

import torch


def distributed_available() -> bool:
    return torch.distributed.is_available() and torch.distributed.is_initialized()


class DistWrapper:
    def __init__(self) -> None:
        self.rank = int(os.environ.get("RANK", 0))
        self.local_rank = int(os.environ.get("LOCAL_RANK", 0))
        self.local_world_size = int(os.environ.get("LOCAL_WORLD_SIZE", 1))
        self.world_size = int(os.environ.get("WORLD_SIZE", 1))
        self.num_nodes = int(self.world_size // self.local_world_size)
        self.node_rank = int(self.rank // self.local_world_size)

    def all_gather_object(self, obj, group=None):
        """Function to gather objects from several distributed processes.
        It is now only used by sync metrics in logger due to security reason.
        """
        if self.world_size > 1 and distributed_available():
            with torch.no_grad():
                obj_list = [None for _ in range(self.world_size)]
                torch.distributed.all_gather_object(obj_list, obj, group=group)
                return obj_list
        else:
            return [obj]


DIST_WRAPPER = DistWrapper()


def traverse_and_aggregate(dict_list, aggregation_func=None):
    """Traverse list of dicts and merge into a single dict with leaf values joined to list."""
    merged_dict = {}
    all_keys = set().union(*dict_list)
    for key in all_keys:
        agg_value = [m[key] for m in dict_list if key in m]

        if isinstance(agg_value[0], dict):
            merged_dict[key] = traverse_and_aggregate(
                agg_value, aggregation_func=aggregation_func
            )
        else:
            if aggregation_func is not None:
                agg_value = aggregation_func(agg_value)
            merged_dict[key] = agg_value

    return merged_dict


def gather_and_merge(metrics, aggregation_func=None):
    """Gather metrics from ddp workers and aggregate leaf metrics."""
    gathered_metrics = DIST_WRAPPER.all_gather_object(metrics)  # list of metrics
    merged_metrics = traverse_and_aggregate(gathered_metrics, aggregation_func)
    return merged_metrics
