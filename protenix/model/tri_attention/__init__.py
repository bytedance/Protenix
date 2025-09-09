# Copyright 2024 ByteDance and/or its affiliates.
# Copyright 2025 Shad Nygren, Virtual Hipster Corporation
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
#
# Modifications by Shad Nygren, Virtual Hipster Corporation (2025):
# - Added fallback implementation for GPUs without Triton support
# - Enables compatibility with consumer GPUs (RTX 3090/4090)
# - Falls back to PyTorch's scaled_dot_product_attention when Triton unavailable

import warnings
import torch

try:
    # Try to import Triton implementation
    import triton
    from protenix.model.tri_attention.op import TriAttention, TriAttentionFunction
    TRITON_AVAILABLE = True
except (ImportError, RuntimeError) as e:
    # Triton not available or not supported on this GPU
    TRITON_AVAILABLE = False
    triton_error = str(e)
    
    # Provide fallback implementation using PyTorch
    class TriAttentionFunction:
        """Fallback implementation using PyTorch's native attention."""
        
        @staticmethod
        def apply(q, k, v, bias1=None, bias2=None, deterministic=False):
            """Apply attention using PyTorch's scaled_dot_product_attention."""
            # q, k, v shape: [B, N, S, H, D]
            B, N, S, H, D = q.shape
            
            # Reshape for attention: [B*N*H, S, D]
            q = q.permute(0, 1, 3, 2, 4).reshape(B * N * H, S, D)
            k = k.permute(0, 1, 3, 2, 4).reshape(B * N * H, S, D)
            v = v.permute(0, 1, 3, 2, 4).reshape(B * N * H, S, D)
            
            # Apply biases if provided
            attn_mask = None
            if bias1 is not None or bias2 is not None:
                # Create attention mask from biases
                # Note: This is a simplified fallback; exact behavior may differ
                attn_mask = 0
                if bias1 is not None:
                    bias1_reshaped = bias1.reshape(B * N * H, 1, S)
                    attn_mask = attn_mask + bias1_reshaped
                if bias2 is not None:
                    bias2_reshaped = bias2.reshape(B * N * H, S, 1)
                    attn_mask = attn_mask + bias2_reshaped
            
            # Use PyTorch's optimized attention
            with torch.backends.cuda.sdp_kernel(
                enable_flash=True,
                enable_math=True,
                enable_mem_efficient=True
            ):
                out = torch.nn.functional.scaled_dot_product_attention(
                    q, k, v,
                    attn_mask=attn_mask,
                    dropout_p=0.0 if deterministic else 0.0,
                    is_causal=False
                )
            
            # Reshape back to original format
            out = out.reshape(B, N, H, S, D).permute(0, 1, 3, 2, 4)
            return out
    
    class TriAttention(torch.nn.Module):
        """Fallback TriAttention module using PyTorch."""
        
        def __init__(self):
            super().__init__()
            warnings.warn(
                f"Triton not available ({triton_error}). Using PyTorch fallback for attention. "
                "This may impact performance but maintains functionality.",
                UserWarning,
                stacklevel=2
            )
        
        def forward(self, q, k, v, bias1=None, bias2=None, deterministic=False):
            return TriAttentionFunction.apply(q, k, v, bias1, bias2, deterministic)

__all__ = ["TriAttention", "TriAttentionFunction", "TRITON_AVAILABLE"]
