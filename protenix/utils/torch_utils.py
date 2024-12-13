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

from contextlib import nullcontext
from typing import Sequence, Union

import numpy as np
import torch
from torch import nn
from torch.nn.parameter import Parameter


def to_device(obj, device):
    """Move tensor or dict of tensors to device"""
    if isinstance(obj, dict):
        for k, v in obj.items():
            if isinstance(v, dict):
                to_device(v, device)
            elif isinstance(v, torch.Tensor):
                obj[k] = obj[k].to(device)
    elif isinstance(obj, torch.Tensor):
        obj = obj.to(device)
    else:
        raise Exception(f"type {type(obj)} not supported")
    return obj


def cdist(a: torch.Tensor, b: torch.Tensor = None):
    # for tensor shape [1, 512 * 14, 3], donot_use_mm_for_euclid_dist mode costs 0.0489s,
    # while use_mm_for_euclid_dist_if_necessary costs 0.0419s on cpu. On GPU there two costs
    # will be neglectible. So there is no need to sacrifice accuracy for speed here.
    return torch.cdist(
        a,
        b if b is not None else a,
        compute_mode="donot_use_mm_for_euclid_dist",
    )


def map_values_to_list(data: dict, recursive: bool = True) -> dict:
    """
    Convert values in a dictionary to lists.

    Args:
        data (dict): The dictionary whose values need to be converted.
        recursive (bool): Whether to recursively convert nested dictionaries. Defaults to True.

    Returns:
        dict: The dictionary with values converted to lists.
    """
    for k, v in data.items():
        if isinstance(v, torch.Tensor):
            if v.dtype == torch.bfloat16:
                v = v.float()
            data[k] = v.cpu().numpy().tolist()
        elif isinstance(v, np.ndarray):
            data[k] = v.tolist()
        elif isinstance(v, dict) and recursive:
            data[k] = map_values_to_list(v, recursive)
    return data


def round_values(data: dict, recursive: bool = True) -> dict:
    """
    Round the values in a dictionary to two decimal places.

    Args:
        data (dict): The dictionary whose values need to be rounded.
        recursive (bool): Whether to recursively round values in nested dictionaries. Defaults to True.

    Returns:
        dict: The dictionary with values rounded to two decimal places.
    """
    for k, v in data.items():
        if isinstance(v, torch.Tensor):
            if v.dtype == torch.bfloat16:
                v = v.float()
            data[k] = np.round(v.cpu().numpy(), 2)
        elif isinstance(v, np.ndarray):
            data[k] = np.round(v, 2)
        elif isinstance(v, list):
            data[k] = list(np.round(np.array(v), 2))
        elif isinstance(v, dict) and recursive:
            data[k] = round_values(v, recursive)
    return data


def autocasting_disable_decorator(disable_casting: bool):
    """
    Decorator to disable autocasting for a function.

    Args:
        disable_casting (bool): If True, disables autocasting; otherwise, uses the default autocasting context.

    Returns:
        function: A decorator that wraps the function with the specified autocasting context.
    """

    def func_wrapper(func):
        def new_func(*args, **kwargs):
            _amp_context = (
                torch.autocast(device_type="cuda", enabled=False)
                if disable_casting
                else nullcontext()
            )
            dtype = torch.float32 if disable_casting else None
            with _amp_context:
                return func(
                    *(
                        v.to(dtype=dtype) if isinstance(v, torch.Tensor) else v
                        for v in args
                    ),
                    **{
                        k: v.to(dtype=dtype) if isinstance(v, torch.Tensor) else v
                        for k, v in kwargs.items()
                    },
                )

        return new_func

    return func_wrapper


def dict_to_tensor(feature_dict: dict) -> dict:
    """
    Convert values in a dictionary to tensors and ensure they have the correct dtype.

    Args:
        feature_dict (dict): The dictionary whose values need to be converted to tensors.

    Returns:
        dict: The dictionary with values converted to tensors and adjusted to the correct dtype.
    """
    for k, v in feature_dict.items():
        if not isinstance(v, torch.Tensor):
            dtype = feature_dict[k].dtype
            feature_dict[k] = torch.tensor(v)

            if dtype in [np.int64, np.int32]:
                feature_dict[k] = feature_dict[k].to(torch.int64)
            elif dtype in [np.float32, np.float64]:
                feature_dict[k] = feature_dict[k].to(torch.float32)

    return feature_dict
