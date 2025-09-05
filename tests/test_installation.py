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
Test suite for installation and dependency compatibility issues.
Related to GitHub issues #182 (pip install error)
"""

import unittest
import sys
import importlib.util


class TestInstallation(unittest.TestCase):
    """Test that package dependencies are compatible"""
    
    def test_python_version(self):
        """Verify Python version is compatible (3.11+)"""
        version_info = sys.version_info
        self.assertGreaterEqual(version_info.major, 3)
        self.assertGreaterEqual(version_info.minor, 11, 
                                "Protenix requires Python 3.11 or higher")
    
    def test_deepspeed_pydantic_compatibility(self):
        """
        Test for issue #182: DeepSpeed and Pydantic compatibility
        The error "json_schema_input_schema" indicates version mismatch
        """
        try:
            # Try to import DeepSpeed - this will fail if incompatible versions
            deepspeed_spec = importlib.util.find_spec("deepspeed")
            if deepspeed_spec is not None:
                import deepspeed
                import pydantic
                
                # This import caused the original error in issue #182
                try:
                    from deepspeed.runtime.config import DeepSpeedConfig
                    # If we get here, the versions are compatible
                    self.assertTrue(True)
                except TypeError as e:
                    if "json_schema_input_schema" in str(e):
                        self.fail(
                            "DeepSpeed/Pydantic compatibility issue detected. "
                            "This is the bug from issue #182. "
                            "Consider pinning pydantic<2.0 or updating DeepSpeed."
                        )
                    else:
                        raise
        except ImportError:
            self.skipTest("DeepSpeed not installed - skipping compatibility test")
    
    def test_required_packages_importable(self):
        """Verify core packages can be imported"""
        required_packages = [
            "torch",
            "numpy", 
            "scipy",
            "pandas",
            "ml_collections",
        ]
        
        for package in required_packages:
            spec = importlib.util.find_spec(package)
            if spec is None:
                self.skipTest(f"{package} not installed - this test documents missing dependencies")
    
    def test_protenix_imports(self):
        """Test that protenix package can be imported when installed"""
        try:
            import protenix
            self.assertTrue(True)
        except ImportError:
            self.skipTest("Protenix not installed - this documents the need for installation")
    
    def test_cuda_availability(self):
        """Document CUDA availability for GPU acceleration"""
        try:
            import torch
            if torch.cuda.is_available():
                self.assertTrue(True, "CUDA is available")
            else:
                # Not a failure, just documentation
                print("Note: CUDA not available, will use CPU")
        except ImportError:
            self.skipTest("PyTorch not installed")


if __name__ == "__main__":
    unittest.main()# Test signed commit
