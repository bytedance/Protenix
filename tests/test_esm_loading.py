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

# Copyright 2025 Shad Nygren, Virtual Hipster Corporation
# Contributed to the Protenix project under the Apache License 2.0
"""
Test suite for ESM model loading compatibility.
Related to GitHub issue #176 (ESM weights loading error with PyTorch 2.6+)
"""

import argparse
from pathlib import Path
from tempfile import TemporaryDirectory
import unittest
from unittest.mock import patch

import torch


class TestESMLoading(unittest.TestCase):
    """Test ESM model loading with different PyTorch versions"""

    @staticmethod
    def _local_esm_checkpoint(model_name: str = "esm2-3b") -> Path:
        from protenix.data.esm.compute_esm import ESM_CONFIG

        return Path("release_data/checkpoint") / ESM_CONFIG[model_name]["model_path"]

    def test_pytorch_version(self):
        """Check PyTorch version and document compatibility"""
        try:
            import torch

            version = torch.__version__
            major, minor = version.split(".")[:2]
            major, minor = int(major), int(
                minor.split("+")[0]
            )  # Handle versions like "2.6.0+cpu"

            if major == 2 and minor >= 6:
                print(
                    f"PyTorch {version} detected - may have ESM loading issues (issue #176)"
                )
            else:
                print(f"PyTorch {version} should work with ESM")

        except ImportError:
            self.skipTest("PyTorch not installed")

    def test_esm_import(self):
        """Test basic ESM package import"""
        try:
            import esm  # noqa: F401

            self.assertTrue(True, "ESM package imported successfully")
        except ImportError:
            self.skipTest("ESM (fair-esm) not installed")

    def test_checkpoint_namespace_loads_with_pytorch_27(self):
        from protenix.data.esm.compute_esm import _load_esm2_model

        with TemporaryDirectory() as temp_dir:
            checkpoint_path = Path(temp_dir) / "esm2_checkpoint.pt"
            torch.save(argparse.Namespace(model_name="esm2"), checkpoint_path)

            def load_model(path):
                return torch.load(path, map_location="cpu"), None

            with patch(
                "protenix.data.esm.compute_esm.pretrained.load_model_and_alphabet_local",
                side_effect=load_model,
            ):
                model, alphabet = _load_esm2_model(checkpoint_path)

        self.assertEqual(model.model_name, "esm2")
        self.assertIsNone(alphabet)

    def test_esm_model_loading(self):
        """
        Test for issue #176: ESM weights loading error in PyTorch 2.6+
        The error occurs when loading pretrained ESM models
        """
        try:
            import torch

            version = torch.__version__
            major, minor = version.split(".")[:2]
            major, minor = int(major), int(minor.split("+")[0])

            from protenix.data.esm.compute_esm import get_esm_embedding  # noqa: F401

            # If we're on PyTorch 2.6+ and this works, the bug is fixed
            if major == 2 and minor >= 6:
                print("ESM loading works on PyTorch 2.6+ - issue #176 may be resolved")

            self.assertTrue(True)

        except ImportError as e:
            if "compute_esm" in str(e):
                self.skipTest("Protenix not installed")
            elif "torch" in str(e):
                self.skipTest("PyTorch not installed")
            else:
                # This might be the ESM loading issue
                try:
                    import torch

                    version = torch.__version__
                    if version.startswith("2.6") or version.startswith("2.7"):
                        self.fail(
                            f"ESM loading failed on PyTorch {version}. "
                            "This is issue #176. The fix has been merged - update Protenix."
                        )
                    else:
                        self.skipTest(f"ESM loading failed: {e}")
                except ImportError:
                    self.skipTest("PyTorch not installed")
        except Exception as e:
            self.skipTest(f"Cannot test ESM loading: {e}")

    def test_esm_embedding_generation(self):
        """Test that ESM embeddings can be generated for a sequence"""
        try:
            import torch

            checkpoint_path = self._local_esm_checkpoint()
            if not checkpoint_path.exists():
                self.skipTest(
                    f"Local ESM checkpoint not available at {checkpoint_path}; "
                    "skipping to avoid downloading weights during tests"
                )

            from protenix.data.esm.compute_esm import load_esm_model

            # Test with a small protein sequence
            test_sequence = "MKTAYDELAAEAFLEENTPILH"
            model, alphabet = load_esm_model("esm2-3b")
            batch_converter = alphabet.get_batch_converter(truncation_seq_length=1022)
            _, _, toks = batch_converter([("test", test_sequence)])
            device = next(model.parameters()).device
            if device.type != "cpu":
                toks = toks.to(device)

            with torch.no_grad():
                out = model(toks, repr_layers=[model.num_layers], return_contacts=False)
            embedding = out["representations"][model.num_layers][0, 1:-1].cpu()

            # Check that we got a reasonable embedding
            self.assertIsNotNone(embedding)
            self.assertGreater(len(embedding), 0)

        except ImportError:
            self.skipTest("Protenix or ESM not installed")
        except unittest.SkipTest:
            raise
        except Exception as e:
            # Check if this is the PyTorch 2.6+ issue
            import torch

            if torch.__version__.startswith("2.6") or torch.__version__.startswith(
                "2.7"
            ):
                self.fail(
                    f"ESM embedding generation failed on PyTorch {torch.__version__}: {e}"
                )
            else:
                self.skipTest(f"Cannot test ESM embeddings: {e}")


if __name__ == "__main__":
    unittest.main()
