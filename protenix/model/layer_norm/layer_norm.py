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

import torch
from torch.nn.parameter import Parameter

sys.path.append(os.path.dirname(__file__))

try:
    fast_layer_norm_cuda_v2 = importlib.import_module("fast_layer_norm_cuda_v2")
except ImportError:
    from protenix.model.layer_norm.torch_ext_compile import compile

    current_dir = os.path.dirname(__file__)
    fast_layer_norm_cuda_v2 = compile(
        name="fast_layer_norm_cuda_v2",
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

        ctx.normalized_shape = normalized_shape
        ctx.eps = eps
        input_ = input.contiguous()

        if weight is None:
            if bias is None:
                output, mean, invvar = fast_layer_norm_cuda_v2.forward_none_affine(
                    input_, ctx.normalized_shape, ctx.eps
                )
            else:
                output, mean, invvar = fast_layer_norm_cuda_v2.forward_with_bias_affine(
                    input_, ctx.normalized_shape, bias.to(d), ctx.eps
                )
        else:
            if bias is None:
                output, mean, invvar = (
                    fast_layer_norm_cuda_v2.forward_with_weight_affine(
                        input_, ctx.normalized_shape, weight.to(d), ctx.eps
                    )
                )
            else:
                output, mean, invvar = fast_layer_norm_cuda_v2.forward_with_both_affine(
                    input_,
                    ctx.normalized_shape,
                    weight.to(d),
                    bias.to(d),
                    ctx.eps,
                )
        ctx.save_for_backward(input_, weight, bias, mean, invvar)
        return output

    @staticmethod
    def backward(ctx, grad_output):
        d = grad_output.dtype
        input_, weight_, bias_, mean, invvar = ctx.saved_tensors
        grad_input = grad_weight = grad_bias = None

        if weight_ is None:
            if bias_ is None:
                grad_input, grad_weight, grad_bias = (
                    fast_layer_norm_cuda_v2.backward_none_affine(
                        grad_output.contiguous(),
                        mean,
                        invvar,
                        input_,
                        ctx.normalized_shape,
                        ctx.eps,
                    )
                )
            else:
                grad_input, grad_weight, grad_bias = (
                    fast_layer_norm_cuda_v2.backward_with_bias_affine(
                        grad_output.contiguous(),
                        mean,
                        invvar,
                        input_,
                        ctx.normalized_shape,
                        bias_.to(dtype=d),
                        ctx.eps,
                    )
                )
        else:
            if bias_ is None:
                grad_input, grad_weight, grad_bias = (
                    fast_layer_norm_cuda_v2.backward_with_weight_affine(
                        grad_output.contiguous(),
                        mean,
                        invvar,
                        input_,
                        ctx.normalized_shape,
                        weight_.to(dtype=d),
                        ctx.eps,
                    )
                )
            else:
                grad_input, grad_weight, grad_bias = (
                    fast_layer_norm_cuda_v2.backward_with_both_affine(
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
        return (
            grad_input,
            None if weight_ is None else grad_weight,
            None if bias_ is None else grad_bias,
            None,
            None,
            None,
        )


class FusedLayerNorm(torch.nn.Module):

    def __init__(
        self, normalized_shape, create_scale=True, create_offset=True, eps=1e-5
    ):
        """
        Args:
            normalized_shape (int or list or torch.Size) input shape from an expected input of size
            create_scale (bool) If set to False, the layer will not learn an additive weight, Default: True
            create_offset (bool) If set to False, the layer will not learn an additive bias, Default: True
            eps (float) a value added to the denominator for numerical stability. Default: 1e-5
        """
        super(FusedLayerNorm, self).__init__()

        if isinstance(normalized_shape, numbers.Integral):
            normalized_shape = (normalized_shape,)
        self.normalized_shape = torch.Size(normalized_shape)
        self.eps = eps
        if create_scale:
            self.weight = Parameter(torch.ones(*normalized_shape))
        else:
            self.weight = None

        if create_offset:
            self.bias = Parameter(torch.zeros(*normalized_shape))
        else:
            self.bias = None

        self.reset_parameters()

    def reset_parameters(self):
        if self.weight is not None:
            torch.nn.init.ones_(self.weight)
        if self.bias is not None:
            torch.nn.init.zeros_(self.bias)

    def forward(self, input):
        return FusedLayerNormAffineFunction.apply(
            input, self.weight, self.bias, self.normalized_shape, self.eps
        )


if __name__ == "__main__":
    dtype = torch.float32
    data = torch.rand(10, 10).cuda().to(dtype=dtype)
    data1 = data * 1
    data.requires_grad = True
    data1.requires_grad = True
    layer_norm = (
        FusedLayerNorm(10, create_scale=True, create_offset=True).cuda().to(dtype=dtype)
    )
    layer_norm_torch = torch.nn.LayerNorm(10).cuda().to(dtype=dtype)
    out = layer_norm(data)
    out1 = layer_norm_torch(data1)
    # print(out - out1)
    loss = out.sum()
    loss.backward()
    loss1 = out1.sum()
    loss1.backward()
    print(data.grad - data1.grad)
    print(layer_norm.weight.grad - layer_norm_torch.weight.grad)
    print(layer_norm.bias.grad - layer_norm_torch.bias.grad)
    print(layer_norm.weight.grad, layer_norm.bias.grad)

    # layer_norm = FusedLayerNorm(10, create_scale=True, create_offset=False).cuda()
    # out = layer_norm(data)
    # loss = out.sum()
    # loss.backward()

    # layer_norm = FusedLayerNorm(10, create_scale=False, create_offset=True).cuda()
    # out = layer_norm(data)
    # loss = out.sum()
    # loss.backward()

    # layer_norm = FusedLayerNorm(10, create_scale=True, create_offset=True).cuda()
    # out = layer_norm(data)
    # loss = out.sum()
    # loss.backward()
