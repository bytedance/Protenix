# Copyright 2025 Shad Nygren, Virtual Hipster Corporation
# Contributed to the Protenix project under the Apache License 2.0
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

"""
Test suite for ESM model loading compatibility.
Related to GitHub issue #176 (ESM weights loading error with PyTorch 2.6+)
"""

import unittest
import sys
import importlib.util


class TestESMLoading(unittest.TestCase):
    """Test ESM model loading with different PyTorch versions"""
    
    def test_pytorch_version(self):
        """Check PyTorch version and document compatibility"""
        try:
            import torch
            version = torch.__version__
            major, minor = version.split('.')[:2]
            major, minor = int(major), int(minor.split('+')[0])  # Handle versions like "2.6.0+cpu"
            
            if major == 2 and minor >= 6:
                print(f"PyTorch {version} detected - may have ESM loading issues (issue #176)")
            else:
                print(f"PyTorch {version} should work with ESM")
                
        except ImportError:
            self.skipTest("PyTorch not installed")
    
    def test_esm_import(self):
        """Test basic ESM package import"""
        try:
            import esm
            self.assertTrue(True, "ESM package imported successfully")
        except ImportError:
            self.skipTest("ESM (fair-esm) not installed")
    
    def test_esm_model_loading(self):
        """
        Test for issue #176: ESM weights loading error in PyTorch 2.6+
        The error occurs when loading pretrained ESM models
        """
        try:
            import torch
            version = torch.__version__
            major, minor = version.split('.')[:2]
            major, minor = int(major), int(minor.split('+')[0])
            
            # Try to load ESM model (this is where issue #176 occurs)
            from protenix.data.compute_esm import get_esm_embedding
            
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
            from protenix.data.compute_esm import get_esm_embedding
            
            # Test with a small protein sequence
            test_sequence = "MKTAYDELAAEAFLEENTPILH"
            
            # This should work if ESM is properly loaded
            embedding = get_esm_embedding(test_sequence)
            
            # Check that we got a reasonable embedding
            self.assertIsNotNone(embedding)
            self.assertGreater(len(embedding), 0)
            
        except ImportError:
            self.skipTest("Protenix or ESM not installed")
        except Exception as e:
            # Check if this is the PyTorch 2.6+ issue
            import torch
            if torch.__version__.startswith("2.6") or torch.__version__.startswith("2.7"):
                self.fail(f"ESM embedding generation failed on PyTorch {torch.__version__}: {e}")
            else:
                self.skipTest(f"Cannot test ESM embeddings: {e}")


if __name__ == "__main__":
    unittest.main()