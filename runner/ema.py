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

from typing import List, Optional

import torch


class EMAWrapper(object):
    """
    A wrapper class for exponential moving average of model weights.

    Args:
        model (torch.nn.Module): The PyTorch model to apply EMA to.
        decay (float): The decay rate for exponential moving average.
        mutable_param_keywords (Optional[List[str]]): List of keywords to identify parameters
            that should be updated by EMA. If None, all parameters are updated.
    """

    def __init__(
        self,
        model: torch.nn.Module,
        decay: float = 0.999,
        mutable_param_keywords: Optional[List[str]] = None,
    ):
        self.model = model
        self.decay = decay
        if mutable_param_keywords is not None:
            self.mutable_param_keywords = [
                s.strip() for s in mutable_param_keywords if s.strip()
            ]
        else:
            self.mutable_param_keywords = None
        self.shadow = {}
        self.backup = {}

    def register(self) -> None:
        """
        Register the model's parameters and create initial shadow copies.
        """
        for name, param in self.model.named_parameters():
            self.shadow[name] = param.data.clone()

    def update(self) -> None:
        """
        Update the shadow copies of the parameters using the current model weights and decay rate.
        """
        for name, param in self.model.named_parameters():
            if self.mutable_param_keywords and not any(
                [keyword in name for keyword in self.mutable_param_keywords]
            ):
                continue
            assert name in self.shadow
            new_average = (1.0 - self.decay) * param.data + self.decay * self.shadow[
                name
            ]
            self.shadow[name] = new_average.clone()

    def apply_shadow(self) -> None:
        """
        Apply the shadow copies (EMA weights) to the model, backing up current weights.
        """
        for name, param in self.model.named_parameters():
            assert name in self.shadow
            self.backup[name] = param.data
            param.data = self.shadow[name]

    def restore(self) -> None:
        """
        Restore the original model weights from the backup.
        """
        for name, param in self.model.named_parameters():
            assert name in self.backup
            param.data = self.backup[name]
        self.backup = {}
