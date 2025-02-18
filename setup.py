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

import sys
from setuptools import find_packages, setup


# Read requirements from the file
with open("requirements.txt") as f:
    install_requires = f.read().splitlines()

# Check if the user specified the CPU option
if '--cpu' in sys.argv:
    # Remove the gpu packages
    try:
        to_drop = [x for x in install_requires if "nvidia" in x or "cuda" in x]
        for x in to_drop:
            install_requires.remove(x)
    except ValueError:
        pass
    # Remove the --cpu option from sys.argv so setuptools doesn't get confused
    sys.argv.remove('--cpu')

setup(
    name="protenix",
    python_requires=">=3.10",
    version="0.4.2",
    description="A trainable PyTorch reproduction of AlphaFold 3.",
    author="Bytedance Inc.",
    url="https://github.com/bytedance/Protenix",
    author_email="ai4s-bio@bytedance.com",
    packages=find_packages(
        exclude=(
            "assets",
            "benchmark",
            "*.egg-info",
        )
    ),
    include_package_data=True,
    package_data={
        "protenix": ["model/layer_norm/kernel/*"],
    },
    install_requires=install_requires,
    license="Apache 2.0 License",
    platforms="manylinux1",
    entry_points={
        "console_scripts": [
            "protenix = runner.batch_inference:protenix_cli",
        ],
    },
)
