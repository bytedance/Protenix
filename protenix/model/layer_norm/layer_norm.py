# Copyright 2024 ByteDance and/or its affiliates.
#
# Copyright 2021- HPC-AI Technology Inc.
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#      http:#www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
import importlib
import numbers
import os
import sys
import time

import torch
from torch.nn.parameter import Parameter

sys.path.append(os.path.dirname(__file__))

try:
    fastfold_layer_norm_cuda = importlib.import_module("fastfold_layer_norm_cuda")
except ImportError:
    from protenix.model.layer_norm.torch_ext_compile import compile

    current_dir = os.path.dirname(__file__)
    fastfold_layer_norm_cuda = compile(
        name="fastfold_layer_norm_cuda",
        sources=[
            os.path.join(f"{current_dir}/kernel", file)
            for file in ["layer_norm_cuda.cpp", "layer_norm_cuda_kernel.cu"]
        ],
        extra_include_paths=[f"{current_dir}/kernel"],
        build_directory=current_dir,
    )


class FusedLayerNormAffineFunction(torch.autograd.Function):

    @staticmethod
    def forward(ctx, input, weight, bias, normalized_shape, eps):
        d = input.dtype
        if d is torch.bfloat16:
            with torch.cuda.amp.autocast(enabled=False):
                ctx.normalized_shape = normalized_shape
                ctx.eps = eps
                input_ = input.contiguous()
                weight_ = weight.contiguous().to(dtype=d)
                bias_ = bias.contiguous().to(dtype=d)
                output, mean, invvar = fastfold_layer_norm_cuda.forward_affine(
                    input_, ctx.normalized_shape, weight_, bias_, ctx.eps
                )
                ctx.save_for_backward(input_, weight_, bias_, mean, invvar)
        else:
            ctx.normalized_shape = normalized_shape
            ctx.eps = eps
            input_ = input.contiguous()
            weight_ = weight.contiguous()
            bias_ = bias.contiguous()
            output, mean, invvar = fastfold_layer_norm_cuda.forward_affine(
                input_, ctx.normalized_shape, weight_, bias_, ctx.eps
            )
            ctx.save_for_backward(input_, weight_, bias_, mean, invvar)

        return output

    @staticmethod
    def backward(ctx, grad_output):
        d = grad_output.dtype
        if d is torch.bfloat16:
            with torch.cuda.amp.autocast(enabled=False):
                input_, weight_, bias_, mean, invvar = ctx.saved_tensors
                grad_input = grad_weight = grad_bias = None
                grad_input, grad_weight, grad_bias = (
                    fastfold_layer_norm_cuda.backward_affine(
                        grad_output.contiguous(),
                        mean,
                        invvar,
                        input_,
                        ctx.normalized_shape,
                        weight_.to(dtype=d),
                        bias_.to(dtype=d),
                        ctx.eps,
                    )
                )
        else:
            input_, weight_, bias_, mean, invvar = ctx.saved_tensors
            grad_input = grad_weight = grad_bias = None
            grad_input, grad_weight, grad_bias = (
                fastfold_layer_norm_cuda.backward_affine(
                    grad_output.contiguous(),
                    mean,
                    invvar,
                    input_,
                    ctx.normalized_shape,
                    weight_,
                    bias_,
                    ctx.eps,
                )
            )

        return grad_input, grad_weight, grad_bias, None, None


class FusedLayerNorm(torch.nn.Module):

    def __init__(self, normalized_shape, eps=1e-5):
        super(FusedLayerNorm, self).__init__()

        if isinstance(normalized_shape, numbers.Integral):
            normalized_shape = (normalized_shape,)
        self.normalized_shape = torch.Size(normalized_shape)
        self.eps = eps
        self.weight = Parameter(torch.ones(*normalized_shape))
        self.bias = Parameter(torch.ones(*normalized_shape))
        self.reset_parameters()

    def reset_parameters(self):
        torch.nn.init.ones_(self.weight)
        torch.nn.init.zeros_(self.bias)

    def forward(self, input):
        return self.kernel_forward(input)

    def kernel_forward(self, input):
        return FusedLayerNormAffineFunction.apply(
            input, self.weight, self.bias, self.normalized_shape, self.eps
        )
