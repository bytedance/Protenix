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

import time
import unittest

import torch
import torch.nn as nn

from protenix.utils.lr_scheduler import AlphaFold3LRScheduler


class SimpleModel(nn.Module):
    def __init__(self):
        super(SimpleModel, self).__init__()
        self.linear = nn.Linear(in_features=1, out_features=1)

    def forward(self, x):
        return self.linear(x)


class TestSchedule(unittest.TestCase):
    def setUp(self):
        self._start_time = time.time()
        return super().setUp()

    def test_af3_lr_schedule(self):
        model = SimpleModel()
        base_lr = 1.8e-3
        optimizer = torch.optim.Adam(
            model.parameters(), lr=base_lr, betas=(0.9, 0.95), eps=1e-8
        )

        scheduler = AlphaFold3LRScheduler(optimizer=optimizer)
        learning_rates = []
        test_steps = 60000
        for step in range(test_steps):
            learning_rates.append(scheduler._get_step_lr(step))
            optimizer.step()
            scheduler.step()
        self.assertEqual(learning_rates[0], 0)
        self.assertEqual(learning_rates[1], 1.8e-6)
        self.assertEqual(learning_rates[1000], 1.8e-3)
        self.assertEqual(learning_rates[50000], 0.95 * 1.8e-3)

    def tearDown(self):
        elapsed_time = time.time() - self._start_time
        print(f"Test {self.id()} took {elapsed_time:.6f}s")


if __name__ == "__main__":
    unittest.main()
