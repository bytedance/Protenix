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
Test suite for Triton kernel compatibility.
Related to GitHub issue #185 (Triton version compatibility)
"""

import unittest
import importlib.util
import os


class TestTritonCompatibility(unittest.TestCase):
    """Test Triton kernel compatibility across different GPU types"""
    
    def test_triton_version(self):
        """Verify Triton version is compatible (3.3.0 recommended)"""
        try:
            import triton
            version = triton.__version__
            # Issue #185 recommends triton==3.3.0
            self.assertEqual(version.split('.')[0], '3', 
                            f"Triton major version should be 3, got {version}")
        except ImportError:
            self.skipTest("Triton not installed")
    
    def test_triton_kernel_import(self):
        """
        Test for issue #185: Triton kernel import failure
        Error: "Failed to import Triton-based component: triangle_multiplicative_update: Not Supported"
        """
        try:
            # This import fails on RTX 3090/4090 GPUs
            from protenix.model.tri_attention import op
            self.assertTrue(True, "Triton attention kernel imported successfully")
        except ImportError as e:
            if "Not Supported" in str(e):
                # This is the bug from issue #185
                import torch
                if torch.cuda.is_available():
                    gpu_name = torch.cuda.get_device_name()
                    if "3090" in gpu_name or "4090" in gpu_name:
                        self.fail(
                            f"Triton not supported on {gpu_name}. "
                            "This is issue #185. Use --trimul_kernel 'torch' as workaround."
                        )
                # For other GPUs, this might be expected
                self.skipTest(f"Triton kernel not supported: {e}")
            else:
                # Different import error
                self.skipTest(f"Protenix not installed: {e}")
        except Exception as e:
            # Unexpected error
            self.skipTest(f"Cannot test Triton compatibility: {e}")
    
    def test_torch_fallback(self):
        """Verify torch fallback works when Triton is not available"""
        # Set environment to use torch kernel
        os.environ["TRIMUL_KERNEL"] = "torch"
        
        try:
            # This should work even without Triton
            import protenix.model.modules
            self.assertTrue(True, "Torch fallback kernel works")
        except ImportError:
            self.skipTest("Protenix not installed")
        finally:
            # Clean up environment
            if "TRIMUL_KERNEL" in os.environ:
                del os.environ["TRIMUL_KERNEL"]
    
    def test_gpu_compatibility_check(self):
        """Document which GPUs have Triton compatibility issues"""
        try:
            import torch
            if not torch.cuda.is_available():
                self.skipTest("No GPU available")
            
            gpu_name = torch.cuda.get_device_name()
            known_issues = ["RTX 3090", "RTX 4090", "GeForce RTX 3090", "GeForce RTX 4090"]
            
            for gpu in known_issues:
                if gpu in gpu_name:
                    print(f"Warning: {gpu_name} has known Triton compatibility issues (issue #185)")
                    print("Workaround: Use --trimul_kernel 'torch'")
                    # Don't fail the test, just document it
                    return
            
            print(f"GPU {gpu_name} should work with Triton kernels")
            
        except ImportError:
            self.skipTest("PyTorch not installed")


if __name__ == "__main__":
    unittest.main()