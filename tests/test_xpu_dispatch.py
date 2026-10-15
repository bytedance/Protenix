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

from types import SimpleNamespace
from unittest import TestCase
from unittest.mock import Mock, patch

import torch

from protenix.config.config import parse_configs
from protenix.utils.device import (
    disable_autocast,
    empty_cache,
    get_available_device_type,
    infer_device_type,
)
from runner.inference import (
    InferenceRunner,
    get_inference_config_defaults,
    run,
    update_gpu_compatible_configs,
)


def inference_configs(triangle_attention=None, triangle_multiplicative=None):
    return SimpleNamespace(
        triangle_attention=triangle_attention,
        triangle_multiplicative=triangle_multiplicative,
        enable_tf32=True,
        dtype="bf16",
    )


class TestXPUDispatch(TestCase):
    def test_cuda_is_preferred_and_keeps_its_kernel_defaults(self):
        with (
            patch("protenix.utils.device.torch.cuda.is_available", return_value=True),
            patch("protenix.utils.device.is_xpu_available", return_value=True),
        ):
            self.assertEqual(get_available_device_type(), "cuda")

        configs = inference_configs()
        with (
            patch("runner.inference.get_available_device_type", return_value="cuda"),
            patch("runner.inference.torch.cuda.is_available", return_value=True),
            patch(
                "runner.inference.torch.cuda.get_device_capability",
                return_value=(8, 0),
            ),
        ):
            update_gpu_compatible_configs(configs)
        self.assertEqual(configs.triangle_attention, "cuequivariance")
        self.assertEqual(configs.triangle_multiplicative, "cuequivariance")
        self.assertTrue(configs.enable_tf32)
        self.assertEqual(configs.dtype, "bf16")

    def test_cpu_keeps_its_kernel_default(self):
        configs = inference_configs()
        with (
            patch("runner.inference.get_available_device_type", return_value="cpu"),
            patch("runner.inference.torch.cuda.is_available", return_value=False),
        ):
            update_gpu_compatible_configs(configs)
        self.assertEqual(configs.triangle_attention, "cuequivariance")
        self.assertEqual(configs.triangle_multiplicative, "cuequivariance")

    def test_xpu_defaults_to_portable_kernels(self):
        configs = inference_configs()
        with patch("runner.inference.get_available_device_type", return_value="xpu"):
            update_gpu_compatible_configs(configs)
        self.assertEqual(configs.triangle_attention, "torch")
        self.assertEqual(configs.triangle_multiplicative, "torch")
        self.assertTrue(configs.enable_tf32)

    def test_xpu_rejects_explicit_non_torch_kernels(self):
        unsupported = {
            "triangle_attention": (
                "cuequivariance",
                "deepspeed",
                "triattention",
            ),
            "triangle_multiplicative": ("cuequivariance",),
        }
        for name, backends in unsupported.items():
            for backend in backends:
                with self.subTest(name=name, backend=backend):
                    configs = inference_configs()
                    setattr(configs, name, backend)
                    with patch(
                        "runner.inference.get_available_device_type",
                        return_value="xpu",
                    ), self.assertRaisesRegex(
                        ValueError, rf"{name}=.*not supported.*use 'torch'"
                    ):
                        update_gpu_compatible_configs(configs)

    def test_xpu_rejects_distributed_inference_before_initialization(self):
        with (
            patch("runner.inference.get_available_device_type", return_value="xpu"),
            patch("runner.inference.DIST_WRAPPER.world_size", 2),
            patch("runner.inference.set_device") as set_device,
            self.assertRaisesRegex(RuntimeError, "single-device inference only"),
        ):
            InferenceRunner(inference_configs())
        set_device.assert_not_called()

    def test_cli_keeps_kernel_provenance_and_selects_xpu_defaults(self):
        with (
            patch("runner.inference.get_available_device_type", return_value="xpu"),
            patch("runner.inference.DIST_WRAPPER.world_size", 1),
            patch(
                "runner.inference.parse_sys_args",
                return_value="--model_name protenix_tiny_default_v0.5.0",
            ),
            patch("runner.inference.download_inference_cache") as download,
            patch("runner.inference.main") as main,
        ):
            run()

        configs = main.call_args.args[0]
        self.assertEqual(configs.triangle_attention, "torch")
        self.assertEqual(configs.triangle_multiplicative, "torch")
        download.assert_called_once()
        self.assertEqual(download.call_args.args[0].triangle_attention, "torch")

    def test_cli_rejects_explicit_cuequivariance_before_download(self):
        with (
            patch("runner.inference.get_available_device_type", return_value="xpu"),
            patch("runner.inference.DIST_WRAPPER.world_size", 1),
            patch(
                "runner.inference.parse_sys_args",
                return_value=(
                    "--model_name protenix_tiny_default_v0.5.0 "
                    "--triangle_attention cuequivariance"
                ),
            ),
            patch("runner.inference.download_inference_cache") as download,
            self.assertRaisesRegex(
                ValueError, "triangle_attention='cuequivariance'.*use 'torch'"
            ),
        ):
            run()
        download.assert_not_called()

    def test_python_entrypoint_uses_xpu_defaults_before_download(self):
        from runner.batch_inference import get_default_runner

        with (
            patch("runner.inference.get_available_device_type", return_value="xpu"),
            patch("runner.inference.DIST_WRAPPER.world_size", 1),
            patch("runner.batch_inference.download_inference_cache") as download,
            patch(
                "runner.batch_inference.InferenceRunner",
                side_effect=lambda configs: configs,
            ),
        ):
            configs = get_default_runner()

        self.assertEqual(configs.triangle_attention, "torch")
        self.assertEqual(configs.triangle_multiplicative, "torch")
        download.assert_called_once()
        self.assertEqual(download.call_args.args[0].triangle_attention, "torch")

    def test_python_entrypoint_rejects_explicit_backend_before_download(self):
        from runner.batch_inference import get_default_runner

        with (
            patch("runner.inference.get_available_device_type", return_value="xpu"),
            patch("runner.inference.DIST_WRAPPER.world_size", 1),
            patch("runner.batch_inference.download_inference_cache") as download,
            self.assertRaisesRegex(
                ValueError, "triangle_attention='cuequivariance'.*use 'torch'"
            ),
        ):
            get_default_runner(triatt_kernel="cuequivariance")
        download.assert_not_called()

    def test_device_inference_prefers_accelerator_over_cpu_metadata(self):
        metadata = {"mask": torch.zeros(1), "features": [torch.device("xpu")]}
        self.assertEqual(infer_device_type(metadata), "xpu")
        self.assertEqual(infer_device_type("xpu:0"), "xpu")

    def test_empty_cache_uses_explicit_device(self):
        with patch("protenix.utils.device.torch.cuda.empty_cache") as cuda_empty_cache:
            empty_cache(torch.device("cuda"))
        cuda_empty_cache.assert_called_once_with()

        xpu_backend = SimpleNamespace(empty_cache=Mock())
        with patch.object(torch, "xpu", xpu_backend, create=True):
            empty_cache(torch.device("xpu"))
        xpu_backend.empty_cache.assert_called_once_with()

        with patch("protenix.utils.device.torch.cuda.empty_cache") as cuda_empty_cache:
            empty_cache(torch.device("cpu"))
        cuda_empty_cache.assert_not_called()

    def test_empty_cache_rejects_non_device_arguments(self):
        with self.assertRaisesRegex(TypeError, "device must be a torch.device"):
            empty_cache(torch.zeros(1))

    def test_disable_autocast_uses_explicit_device(self):
        with patch("protenix.utils.device.torch.amp.autocast") as autocast:
            disable_autocast(torch.device("xpu"))
        autocast.assert_called_once_with("xpu", enabled=False)

        with patch("protenix.utils.device.torch.amp.autocast") as autocast:
            disable_autocast(torch.device("cuda"))
        autocast.assert_called_once_with("cuda", enabled=False)

        with patch("protenix.utils.device.torch.amp.autocast") as autocast:
            disable_autocast(torch.device("cpu"))
        autocast.assert_not_called()

    def test_disable_autocast_rejects_non_device_arguments(self):
        with self.assertRaisesRegex(TypeError, "device must be a torch.device"):
            disable_autocast(torch.zeros(1))

    def test_inference_defaults_mark_kernel_choices_as_automatic(self):
        configs = parse_configs(
            get_inference_config_defaults(),
            arg_str="--model_name protenix_tiny_default_v0.5.0",
            fill_required_with_null=True,
        )
        self.assertIsNone(configs.triangle_attention)
        self.assertIsNone(configs.triangle_multiplicative)

        explicit = parse_configs(
            get_inference_config_defaults(),
            arg_str=(
                "--model_name protenix_tiny_default_v0.5.0 "
                "--triangle_attention cuequivariance"
            ),
            fill_required_with_null=True,
        )
        self.assertEqual(explicit.triangle_attention, "cuequivariance")
        self.assertIsNone(explicit.triangle_multiplicative)
