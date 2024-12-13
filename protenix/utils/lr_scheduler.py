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

import math
import warnings

import torch
from torch.optim.lr_scheduler import LRScheduler


class CosineAnnealingWithWarmup(LRScheduler):
    def __init__(
        self,
        optimizer: torch.optim.Optimizer,
        warmup_steps: int,
        decay_steps: int,
        lr: float,
        min_lr: float,
        last_epoch: int = -1,
        verbose: bool = False,
    ):
        self.warmup_steps = warmup_steps
        self.decay_steps = decay_steps
        self.lr = lr
        self.min_lr = min_lr
        super().__init__(optimizer, last_epoch, verbose)

    def _get_step_lr(self, step):
        if step <= self.warmup_steps:
            return (step + 1) / (self.warmup_steps + 1) * self.lr
        elif step >= self.decay_steps:
            return self.min_lr
        else:
            decay_ratio = (step - self.warmup_steps) / (
                self.decay_steps - self.warmup_steps
            )
            assert 0 <= decay_ratio <= 1
            coff = 0.5 * (1.0 + math.cos(math.pi * decay_ratio))
            return self.min_lr + coff * (self.lr - self.min_lr)

    def get_lr(self):
        if not self._get_lr_called_within_step:
            warnings.warn(
                "To get the last learning rate computed by the scheduler, "
                "please use `get_last_lr()`.",
                UserWarning,
            )
        return [
            self._get_step_lr(self.last_epoch) for group in self.optimizer.param_groups
        ]

    def _get_closed_form_lr(self):
        return [self._get_step_lr(self.last_epoch) for base_lr in self.base_lrs]


# The Alphafold3 Learning Rate Scheduler As in 5.4
class AlphaFold3LRScheduler(LRScheduler):
    def __init__(
        self,
        optimizer: torch.optim.Optimizer,
        last_epoch: int = -1,
        verbose: bool = False,
        warmup_steps: int = 1000,
        lr: float = 1.8e-3,
        decay_every_n_steps: int = 50000,
        decay_factor: float = 0.95,
    ) -> None:
        self.warmup_steps = warmup_steps
        self.decay_steps = decay_every_n_steps
        self.lr = lr
        self.decay_factor = decay_factor
        super(AlphaFold3LRScheduler, self).__init__(
            optimizer=optimizer, last_epoch=last_epoch, verbose=verbose
        )

    def _get_step_lr(self, step):
        if step <= self.warmup_steps:
            lr = step / self.warmup_steps * self.lr
        else:
            decay_count = step // self.decay_steps
            lr = self.lr * (self.decay_factor**decay_count)
        return lr

    def get_lr(self) -> list[float]:
        if not self._get_lr_called_within_step:
            warnings.warn(
                "To get the last learning rate computed by the scheduler, "
                "please use `get_last_lr()`.",
                UserWarning,
            )
        return [
            self._get_step_lr(self.last_epoch) for group in self.optimizer.param_groups
        ]


def get_lr_scheduler(
    configs, optimizer: torch.optim.Optimizer, **kwargs
) -> torch.optim.lr_scheduler.LRScheduler:
    """
    Get the learning rate scheduler based on the configuration.

    Args:
        configs: Configuration object containing scheduler settings.
        optimizer (torch.optim.Optimizer): The optimizer to which the scheduler will be attached.
        **kwargs: Additional keyword arguments to be passed to the scheduler.

    Returns:
        torch.optim.lr_scheduler.LRScheduler: The learning rate scheduler.

    Raises:
        ValueError: If the specified learning rate scheduler is invalid.
    """
    if configs.lr_scheduler == "af3":
        lr_scheduler = AlphaFold3LRScheduler(
            optimizer, **configs.af3_lr_scheduler, **kwargs
        )
    elif configs.lr_scheduler == "cosine_annealing":
        lr_scheduler = CosineAnnealingWithWarmup(
            optimizer,
            configs.warmup_steps,
            configs.max_steps,
            configs.lr,
            configs.lr * configs.min_lr_ratio,
            **kwargs,
        )
    elif configs.lr_scheduler == "constant":
        lr_scheduler = torch.optim.lr_scheduler.ConstantLR(
            optimizer,
            factor=1.0,
            total_iters=configs.max_steps,
            **kwargs,
        )
    else:
        raise ValueError(f"Invalid lr scheduler: [{configs.lr_scheduler}]")
    return lr_scheduler
