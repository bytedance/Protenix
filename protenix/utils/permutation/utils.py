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
import time

import torch


class Checker:

    @staticmethod
    def is_permutation(x: torch.Tensor):
        """
        Checks if the input tensor `x` is a permutation of integers from 0 to N-1.

        Args:
            x (torch.Tensor): A 1D tensor of size [N].
        """
        assert x.dim() == 1
        N = x.size(0)
        assert torch.equal(torch.sort(x)[0], torch.arange(N, device=x.device))

    @staticmethod
    def are_permutations(x: torch.Tensor, dim=-1):
        """
        Checks if slices along the specified dimension in `x` are permutations of integers from 0 to N-1.

        Args:
            x (torch.Tensor): A tensor with any number of dimensions, containing slices of size N along `dim`.
            dim (int, optional): The dimension along which to check for permutations. Defaults to -1.
        """
        assert x.dim() > 0

        N = x.size(dim)
        # Create a view of x that moves the specified dimension to -1
        x = x.transpose(dim, -1).contiguous()
        x = x.reshape(-1, N)
        expected = torch.arange(N, device=x.device)
        for i in range(x.size(0)):
            Checker.is_permutation(x[i])

    @staticmethod
    def contains_identity(x: torch.Tensor, dim=-1):
        """
        Check if x contains the identity permutation

        Args:
            x (torch.Tensor): A tensor with any number of dimensions, containing slices of size N along `dim`.
            dim (int, optional): The dimension along which to check for permutations. Defaults to -1.
        """
        assert x.dim() > 0

        N = x.size(dim)
        # Create a view of x that moves the specified dimension to -1
        x = x.transpose(dim, -1).contiguous()
        x = x.reshape(-1, N)
        expected = torch.arange(N, device=x.device).unsqueeze(dim=0)
        assert (x == expected).all(dim=-1).any()

    @staticmethod
    def not_contain_identity(x: torch.Tensor, dim=-1):
        """
        Check if x does not contain the identity permutation

        Args:
            x (torch.Tensor): A tensor with any number of dimensions, containing slices of size N along `dim`.
            dim (int, optional): The dimension along which to check for permutations. Defaults to -1.
        """
        assert x.dim() > 0

        N = x.size(dim)
        # Create a view of x that moves the specified dimension to -1
        x = x.transpose(dim, -1).contiguous()
        x = x.reshape(-1, N)
        expected = torch.arange(N, device=x.device).unsqueeze(dim=0)
        assert not (x == expected).all(dim=-1).any()

    @staticmethod
    def batch_permute(perm: torch.Tensor, x: torch.Tensor, x_permuted: torch.Tensor):
        """
        Args:
            perm (torch.Tensor):
                [..., N]
            x (torch.Tensor):
                [N, batch_dims_x]
            x_permuted (torch.Tensor):
                [..., N, batch_dims_x]
        """
        batch_shape = perm.shape[:-1]
        N = perm.size(-1)
        assert x.size(0) == N
        perm = perm.view(-1, N)
        permuted_x = [x[perm[i]] for i in range(len(perm))]
        permuted_x = torch.stack(permuted_x, dim=0)  # [-1, N, batch_dims_x]
        target_shape = batch_shape + (N,) + x.shape[1:]
        assert torch.allclose(permuted_x.reshape(target_shape), x_permuted)


def save_permutation_error(data, error_dir: str = None, max_cases: int = 50):
    """
    Saves the permutation error data to a specified directory.

    Args:
        data: The data to be saved.
        error_dir (str): The directory where the error data should be saved.
        max_cases (int): The maximum number of error cases to save.

    Raises:
        Exception: If an error occurs while saving the data, the exception is caught and printed.
    """
    if error_dir is None:
        return

    # error_dir = os.path.join(self.error_dir, dir_name)
    os.makedirs(error_dir, exist_ok=True)

    if len(os.listdir(error_dir)) >= max_cases:
        # Only record the first {max_cases} error cases for debug
        return

    filename = "T_" + time.strftime("%Y%m%d_%H%M%S") + ".pt"
    fpath = os.path.join(error_dir, filename)
    if not os.path.exists(fpath):
        try:
            torch.save(data, fpath)
        except Exception as e:
            print(f"Exception occurrs in save_permutation_error: {e}")
