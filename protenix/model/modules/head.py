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

import torch
import torch.nn as nn
from torch.nn import Linear


# Adapted From openfold.model.heads
class DistogramHead(nn.Module):
    """Implements Algorithm 1 [Line17] in AF3
    Computes a distogram probability distribution.
    For use in computation of distogram loss, subsection 1.9.8 (AF2)
    """

    def __init__(self, c_z: int = 128, no_bins: int = 64) -> None:
        """
        Args:
            c_z (int, optional): hidden dim [for pair embedding]. Defaults to 128.
            no_bins (int, optional): Number of distogram bins. Defaults to 64.
        """
        super(DistogramHead, self).__init__()

        self.c_z = c_z
        self.no_bins = no_bins

        self.linear = Linear(in_features=self.c_z, out_features=self.no_bins)

    def forward(self, z: torch.Tensor) -> torch.Tensor:  # [*, N, N, C_z]
        """
        Args:
            z (torch.Tensor): pair embedding
                [*, N_token, N_token, C_z]

        Returns:
            torch.Tensor: distogram probability distribution
                [*, N_token, N_token, no_bins]
        """
        # [*, N, N, no_bins]
        logits = self.linear(z)
        logits = logits + logits.transpose(-2, -3)
        return logits
