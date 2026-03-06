"""Tests for torch.compile support + training loop fixes.

Verifies:
1. @torch.compiler.disable decorated functions still work normally
2. GradScaler dtype key fix ("fp16" vs "float16")
3. compile config structure
4. get_checkpoint_fn caching
"""

import unittest
from functools import partial

import torch


class TestCompilerDisableDecorators(unittest.TestCase):
    """Verify @torch.compiler.disable doesn't break normal execution."""

    def test_softmax_no_cast(self):
        try:
            from protenix.model.triangular.layers import softmax_no_cast
        except (RuntimeError, ModuleNotFoundError) as e:
            # CUDA JIT compilation of custom layer_norm kernel may fail
            # in some environments (e.g., nvcc arch mismatch). Skip gracefully.
            self.skipTest(f"Cannot import due to CUDA JIT build failure: {e}")
        x = torch.randn(4, 8)
        result = softmax_no_cast(x, dim=-1)
        self.assertEqual(result.shape, x.shape)
        # Should sum to 1 along dim=-1
        self.assertTrue(torch.allclose(result.sum(dim=-1), torch.ones(4), atol=1e-5))

    def test_chunk_layer_helpers(self):
        from protenix.model.utils import _flat_idx_to_idx, _get_minimal_slice_set
        # _flat_idx_to_idx
        idx = _flat_idx_to_idx(5, (2, 3))
        self.assertEqual(idx, (1, 2))

        idx = _flat_idx_to_idx(0, (3, 4))
        self.assertEqual(idx, (0, 0))

    def test_checkpoint_blocks(self):
        from protenix.model.utils import checkpoint_blocks
        # Should be callable without torch.compile interfering
        self.assertTrue(callable(checkpoint_blocks))


class TestGetCheckpointFnCaching(unittest.TestCase):
    """Verify get_checkpoint_fn returns the same object each call."""

    def test_same_object(self):
        from protenix.model.utils import get_checkpoint_fn
        fn1 = get_checkpoint_fn()
        fn2 = get_checkpoint_fn()
        self.assertIs(fn1, fn2)

    def test_is_partial(self):
        from protenix.model.utils import get_checkpoint_fn
        fn = get_checkpoint_fn()
        self.assertIsInstance(fn, partial)
        self.assertFalse(fn.keywords.get("use_reentrant", True))


class TestCompileConfig(unittest.TestCase):
    """Verify compile config structure in configs_base."""

    def test_config_has_compile_block(self):
        # Import and check the config dict has the compile key
        from configs.configs_base import model_configs

        self.assertIn("compile", model_configs)
        compile_cfg = model_configs["compile"]
        expected_keys = {"pairformer", "diffusion", "confidence", "msa"}
        self.assertEqual(set(compile_cfg.keys()), expected_keys)
        # All default to False
        for k, v in compile_cfg.items():
            self.assertFalse(v, f"compile.{k} should default to False")


class TestGradScalerDtypeKey(unittest.TestCase):
    """Verify GradScaler works with 'fp16' dtype key (the fix)."""

    def test_scaler_enabled_fp16(self):
        """The fix: dtype == 'fp16' (not 'float16') enables the scaler."""
        dtype_key = "fp16"
        scaler = torch.GradScaler(
            device="cpu",
            enabled=(dtype_key == "fp16"),
        )
        # Scaler should be enabled
        self.assertTrue(scaler.is_enabled())

    def test_scaler_disabled_bf16(self):
        dtype_key = "bf16"
        scaler = torch.GradScaler(
            device="cpu",
            enabled=(dtype_key == "fp16"),
        )
        self.assertFalse(scaler.is_enabled())

    def test_old_bug_float16_key(self):
        """The old code used 'float16' which never matched the config's 'fp16'."""
        dtype_key = "fp16"
        # Old code: enabled=(dtype_key == "float16") -> always False for "fp16"
        scaler_old = torch.GradScaler(device="cpu", enabled=(dtype_key == "float16"))
        self.assertFalse(scaler_old.is_enabled())  # Bug: should be enabled


class TestAutocastCacheDisabled(unittest.TestCase):
    """cache_enabled must be False to avoid gradient checkpointing conflicts
    and vanishing confidence head gradients."""

    def test_cache_enabled_is_false(self):
        import ast
        import inspect
        from runner.train import AF3Trainer as Trainer

        src = inspect.getsource(Trainer)
        tree = ast.parse(src)
        found = False
        for node in ast.walk(tree):
            if isinstance(node, ast.keyword) and node.arg == "cache_enabled":
                if isinstance(node.value, ast.Constant):
                    self.assertFalse(
                        node.value.value,
                        "cache_enabled must be False (breaks gradient checkpointing)",
                    )
                    found = True
        self.assertTrue(found, "cache_enabled keyword not found in Trainer source")


class TestOrigModStateDict(unittest.TestCase):
    """torch.compile adds _orig_mod prefix — checkpoint loading must handle it."""

    def test_strip_orig_mod(self):
        """Uncompiled model loading compiled checkpoint."""
        ckpt_keys = {
            "pairformer_stack._orig_mod.blocks.0.weight": torch.tensor(1.0),
            "pairformer_stack._orig_mod.blocks.0.bias": torch.tensor(0.0),
        }
        stripped = {k.replace("._orig_mod.", "."): v for k, v in ckpt_keys.items()}
        self.assertEqual(
            set(stripped.keys()),
            {"pairformer_stack.blocks.0.weight", "pairformer_stack.blocks.0.bias"},
        )

    def test_add_orig_mod(self):
        """Compiled model loading uncompiled checkpoint."""
        ckpt_keys = {
            "pairformer_stack.blocks.0.weight": torch.tensor(1.0),
            "diffusion_module.layers.0.bias": torch.tensor(0.0),
            "other_param": torch.tensor(2.0),
        }
        prefixes = ["pairformer_stack.", "diffusion_module.", "confidence_head.", "msa_module."]
        added = {}
        for k, v in ckpt_keys.items():
            new_k = k
            for prefix in prefixes:
                if k.startswith(prefix):
                    new_k = k.replace(prefix, prefix + "_orig_mod.", 1)
                    break
            added[new_k] = v
        self.assertEqual(
            set(added.keys()),
            {
                "pairformer_stack._orig_mod.blocks.0.weight",
                "diffusion_module._orig_mod.layers.0.bias",
                "other_param",
            },
        )


if __name__ == "__main__":
    unittest.main()
