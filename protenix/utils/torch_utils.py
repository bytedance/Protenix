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


def grad_norm(params):
    total_norm = 0.0
    for p in params:
        if p.grad is None:
            continue
        param_norm = p.grad.data.norm(2)
        total_norm += param_norm.item() ** 2
    total_norm = total_norm ** (1.0 / 2)
    return total_norm


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


def detach_if(t: torch.Tensor, detach: bool):
    if detach:
        return t.detach()
    return t


def cdist(a: torch.Tensor, b: torch.Tensor = None):
    # for tensor shape [1, 512 * 14, 3], donot_use_mm_for_euclid_dist mode costs 0.0489s,
    # while use_mm_for_euclid_dist_if_necessary costs 0.0419s on cpu. On GPU there two costs
    # will be neglectible. So there is no need to sacrifice accuracy for speed here.
    return torch.cdist(
        a,
        b if b is not None else a,
        compute_mode="donot_use_mm_for_euclid_dist",
    )


def batch_avg_with_mask(
    value: torch.Tensor,
    mask: torch.Tensor,
    avg_dim: Union[int, tuple[int]] = None,
    batch_reduction: str = "mean",
    eps: float = 1e-12,
):
    """Average values with mask.
    Args:
        value: tensor of shape [BS, ...]
        mask: tensor with same shape and type of value, 1 means valid, 0 means maksed
        avg_dim: dimensions to apply average, if None, all dims excluding BS dim will be averaged
        batch_reduction: mean/sum/none, reduction operation applied on BS dim
    """
    if avg_dim is None:
        avg_dim = tuple(range(1, len(value.shape)))
    avg = (value * mask).sum(dim=avg_dim) / (mask.sum(dim=avg_dim) + eps)
    if batch_reduction == "mean":
        return avg.mean()
    elif batch_reduction == "sum":
        return avg.sum()
    elif batch_reduction == "none":
        return avg
    else:
        raise Exception(f"Invalid batch_reduction: {batch_reduction}")


def eye_mask(L, device=None, opposite=False):
    if opposite:
        return 1.0 - torch.eye(L, device=device)
    else:
        torch.eye(L, device=device)


def glorot_uniform(t):
    if len(t.size()) == 2:
        fan_in, fan_out = t.size()
    elif len(t.size()) == 3:
        # out_ch, in_ch, kernel for Conv 1
        fan_in = t.size()[1] * t.size()[2]
        fan_out = t.size()[0] * t.size()[2]
    else:
        fan_in = np.prod(t.size())
        fan_out = np.prod(t.size())

    limit = np.sqrt(6.0 / (fan_in + fan_out))
    t.uniform_(-limit, limit)


def _param_init(m, bias="zero"):
    if isinstance(m, Parameter):
        glorot_uniform(m.data)
    elif isinstance(m, nn.Linear):
        if m.bias is not None:
            if bias == "zero":
                m.bias.data.zero_()
            else:
                assert bias == "normal"
                m.bias.data.normal_()
        glorot_uniform(m.weight.data)


def weights_init(m, bias="zero"):
    for p in m.modules():
        if isinstance(p, nn.ParameterList) or isinstance(p, nn.ModuleList):
            for pp in p:
                _param_init(pp, bias)
        else:
            _param_init(p, bias)

    for name, p in m.named_parameters():
        if "." not in name:  # top-level parameters
            _param_init(p, bias)


def permute_last_dims(t: torch.Tensor, dims: Sequence[int]):
    """Permute tensor on last dims, all other dims are kept unchanged.

    Args:
        t (torch.Tensor): Input tensor with at least len(dims) dimensions.
        dims: The desired ordering of dimensions, here all values should be < 0, i.e. (-1, -2) means permute last two dims.
    """
    num_dims = len(t.shape)
    prefix_dims = list(range(num_dims - len(dims)))
    last_dims = [num_dims + d for d in dims]
    return torch.permute(t, prefix_dims + last_dims)


def flatten_tensors(tensors) -> torch.Tensor:
    """Flatten a list of tensors into a single 1D tensor."""
    return torch.cat([t.view(-1) for t in tensors], dim=0)


def unflatten_tensors(flat_tensor, shapes):
    """Unflatten a 1D tensor into a list of tensors with given shapes."""
    tensors = []
    offset = 0
    for shape in shapes:
        numel = shape.numel()
        tensors.append(flat_tensor[offset : offset + numel].view(shape))
        offset += numel
    return tensors


def map_values_to_list(data, recursive=True):
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


def round_values(data, recursive=True):
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


def autocasting_disable_decorator(disable_casting):
    def func_wrapper(func):
        def new_func(*args, **kwargs):
            _amp_context = (
                torch.autocast(device_type="cuda", enabled=False)
                if disable_casting
                else nullcontext()
            )

            # Helper function to conditionally cast tensors
            def conditioned_cast(tensor):
                if (
                    disable_casting
                    and isinstance(tensor, torch.Tensor)
                    and torch.is_floating_point(tensor)
                ):
                    return tensor.to(dtype=torch.float32)
                return tensor

            with _amp_context:
                return func(
                    *(conditioned_cast(v) for v in args),
                    **{k: conditioned_cast(v) for k, v in kwargs.items()},
                )

        return new_func

    return func_wrapper


def dict_to_tensor(feature_dict):
    for k, v in feature_dict.items():
        if not isinstance(v, torch.Tensor):
            dtype = feature_dict[k].dtype
            feature_dict[k] = torch.tensor(v)

            if dtype in [np.int64, np.int32]:
                feature_dict[k] = feature_dict[k].to(torch.int64)
            elif dtype in [np.float32, np.float64]:
                feature_dict[k] = feature_dict[k].to(torch.float32)

    return feature_dict


def collate_fn_identity(x):
    return x


def collate_fn_first(x):
    return x[0]
