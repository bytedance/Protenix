import torch
import unittest
import numpy as np
from typing import Tuple
import gc
import os
import warnings

from protenix.openfold_local.model.triangular_multiplicative_update import (
    TriangleMultiplicationOutgoing,
    TriangleMultiplicationIncoming,
)


class TestTriangularMultiplicativeUpdateMemory(unittest.TestCase):
    """Test cases comparing memory consumption of different implementations"""
    
    def setUp(self):
        self.device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
        # Test configurations
        self.train_test_configs = [
            # (batch_size, seq_len, hidden_dim)
            (1, 128, 128),
            (1, 512, 128),
            (1, 768, 128),
            (1, 1024, 128),
        ]
        
        self.eval_test_configs = [
            # (batch_size, seq_len, hidden_dim)
            (1, 128, 128),
            (1, 512, 128),
            (1, 768, 128),
            (1, 1024, 128),
            (1, 2048, 128),
            (1, 4096, 128),
        ]
        self.directions = ['outgoing', 'incoming']
        # Reduced iterations for expensive operations like backward pass
        self.num_warmup = 5
        self.num_iterations = 20
        self.atol = 1e-5
        self.rtol = 1e-5
        self.precision = torch.bfloat16
        
        self.ln_init_incoming = torch.nn.LayerNorm(128)
        self.ln_init_incoming.weight.data = torch.tensor(
            [1.1069, 0.1331, 0.8551, 0.9301, 1.3930, 1.2957, 1.1911, 1.2627, 1.0659,
            0.8516, 1.4179, 1.2899, 1.1479, 0.5468, 1.1895, 1.5419, 1.4297, 1.2888,
            1.3909, 0.2827, 1.2673, 1.2662, 0.9909, 0.7268, 1.0813, 1.2905, 1.2443,
            1.1072, 1.0935, 0.8666, 1.4247, 1.5937, 0.9327, 1.4187, 1.0947, 1.3388,
            0.7701, 0.9646, 1.5998, 0.8428, 1.0077, 0.8284, 0.8372, 1.2660, 1.2170,
            1.1748, 1.9523, 1.4707, 0.8892, 0.9333, 1.3764, 1.1505, 1.3801, 1.1086,
            1.5976, 0.8450, 1.1612, 0.9779, 0.8710, 0.8558, 1.4008, 0.9696, 0.9381,
            1.1749, 1.3342, 0.7880, 0.7195, 1.3028, 1.5130, 1.1993, 1.0383, 1.0001,
            1.0637, 1.1721, 0.8464, 0.6833, 1.2382, 1.3218, 1.1659, 1.0175, 1.0792,
            0.9469, 1.0529, 0.9846, 1.0648, 1.5787, 1.4551, 0.9300, 1.2841, 1.2093,
            1.2799, 1.4551, 1.4075, 1.5533, 1.6966, 1.3683, 1.1345, 0.8771, 1.5656,
            1.2981, 0.5383, 1.6160, 1.0285, 0.9068, 1.4749, 1.1919, 1.0276, 1.2010,
            1.3667, 1.3054, 1.2131, 1.2671, 0.8501, 0.9664, 1.1108, 1.0416, 0.9263,
            1.0818, 1.0534, 1.1937, 0.8772, 1.0710, 1.1215, 1.2230, 1.3231, 1.1420,
            1.0391, 1.5955])
        self.ln_init_incoming.bias.data = torch.tensor(
            [ 0.0331,  1.1233, -0.2818,  0.1086, -0.0717,  0.2010, -0.1568,  0.0275,
            0.0991, -0.3981,  0.1851, -0.0672,  0.0486,  0.5538, -0.2101, -0.3011,
            0.0057,  0.1074, -0.0488, -0.8954, -0.2224,  0.2573, -0.1943, -0.6369,
            0.0229, -0.1514,  0.0476,  0.2517,  0.1683, -0.2041, -0.2544, -0.0727,
            0.3481,  0.0924,  0.0208,  0.1356,  0.3388,  0.1527,  0.0581,  0.0080,
            0.3269, -0.3701, -0.0627, -0.1449, -0.0015, -0.0848,  0.1265,  0.2727,
            -0.0854, -0.3870, -0.2771, -0.1739, -0.0833,  0.1799,  0.1896, -0.3888,
            0.1363, -0.2650, -0.3566, -0.0974, -0.3219,  0.1624, -0.2751,  0.0530,
            -0.2622,  0.0575, -0.2851, -0.1549,  0.0772, -0.1936,  0.1195,  0.3170,
            -0.0882,  0.1159, -0.5661,  0.1842, -0.3676,  0.1350, -0.4965,  0.2311,
            0.0482, -0.0558,  0.0491, -0.2848, -0.1151,  0.1570, -0.3073, -0.1670,
            0.0454,  0.1157,  0.0581, -0.2814, -0.1413,  0.0293, -0.0091, -0.1387,
            0.0845, -0.6806, -0.0413,  0.2504,  0.7005,  0.0624,  0.0920, -0.3810,
            0.2874, -0.2491, -0.0779, -0.3682,  0.0235, -0.1377, -0.2093, -0.3319,
            -0.0830,  0.2899,  0.2727,  0.2091, -0.3262, -0.1919,  0.2730, -0.0959,
            0.0283,  0.1232,  0.0608, -0.2017,  0.1202,  0.2317, -0.2799,  0.0505])
        
        self.ln_init_outgoing = torch.nn.LayerNorm(128)
        self.ln_init_outgoing.weight.data = torch.tensor(
            [1.0365, 0.1086, 0.9429, 1.2368, 0.9882, 0.7943, 0.8879, 1.1810, 1.3339,
            0.9879, 1.1907, 1.4831, 1.2308, 0.4383, 1.3593, 1.2732, 1.3220, 1.1092,
            1.0285, 0.2099, 1.4191, 1.3977, 1.4036, 0.7775, 1.4555, 1.4130, 1.1095,
            0.8518, 1.3124, 1.2264, 1.1767, 1.3703, 0.8941, 1.6451, 1.3698, 1.0190,
            0.9669, 1.2591, 1.5631, 1.4062, 1.1375, 1.3069, 0.8156, 1.3431, 1.2776,
            0.8445, 1.0338, 0.9814, 1.2690, 0.6945, 0.6919, 1.1400, 1.3023, 1.0900,
            1.2318, 0.8767, 1.3103, 1.0040, 1.0784, 1.0948, 1.3041, 1.1809, 1.0819,
            1.3783, 1.4137, 0.9668, 0.7806, 1.0276, 1.3490, 1.3316, 1.0288, 1.0727,
            1.2222, 0.9263, 0.6448, 1.1502, 0.7552, 1.3718, 1.4600, 0.7817, 1.3911,
            0.9506, 1.4980, 0.9404, 1.4928, 1.5902, 1.3479, 1.3230, 0.9561, 1.0794,
            1.3996, 1.0368, 1.3119, 1.2974, 1.3994, 1.3613, 1.0881, 0.8575, 1.2145,
            0.6367, 0.4945, 1.8550, 1.3011, 0.8541, 1.5343, 1.3083, 1.1008, 0.8480,
            1.4317, 0.8955, 1.2785, 1.1039, 1.2231, 1.1548, 0.8755, 1.3220, 0.4856,
            1.4148, 1.0072, 1.3423, 1.0607, 0.8498, 1.5235, 1.4233, 1.3585, 0.7820,
            1.3403, 1.2852])
        self.ln_init_outgoing.bias.data = torch.tensor(
            [ 0.0074,  0.9458, -0.3074,  0.1211, -0.0117,  0.0597, -0.2677,  0.1065,
            0.2065, -0.3813,  0.0619, -0.2894, -0.0258,  0.4084, -0.2881,  0.0685,
            0.0665,  0.1907, -0.1177, -0.8015, -0.2551, -0.0175, -0.2045, -0.5616,
            -0.2281, -0.2304,  0.1300,  0.0422,  0.2689, -0.1146, -0.2879,  0.0699,
            0.2924,  0.0269, -0.0904,  0.0585,  0.3193,  0.2071, -0.0303,  0.0751,
            0.3909, -0.3810, -0.0117, -0.2193,  0.1019,  0.0602, -0.0359, -0.0950,
            -0.0374, -0.2458, -0.0839, -0.0484, -0.0103,  0.0117, -0.0335, -0.4499,
            0.0925, -0.2208, -0.3235, -0.0900, -0.2144,  0.2907, -0.2043, -0.0276,
            -0.0392, -0.0309, -0.4329, -0.0664,  0.0637, -0.0523,  0.0941,  0.3192,
            -0.1626,  0.0090, -0.5427,  0.1415, -0.2097,  0.0525, -0.2909,  0.3014,
            0.0924,  0.0011,  0.1771, -0.2014, -0.1694,  0.0220, -0.1532, -0.1851,
            0.0972,  0.0388,  0.0920, -0.2167, -0.1332, -0.0824, -0.1111, -0.0307,
            0.1437, -0.6457,  0.0145,  0.1429,  0.5746,  0.0355,  0.1116, -0.3221,
            0.1196, -0.2788, -0.1722, -0.1533,  0.0254,  0.0170, -0.0634, -0.3098,
            -0.1110,  0.3204,  0.1493,  0.3524, -0.3108, -0.3455,  0.2530, -0.0674,
            0.2484,  0.1975,  0.0975, -0.1844,  0.1865,  0.1264, -0.2840, -0.0640])
        
    def _clear_memory(self):
        """Clear GPU memory cache"""
        gc.collect()
        if self.device.type == 'cuda':
            torch.cuda.empty_cache()
            torch.cuda.synchronize()
    
    def _get_memory_usage(self) -> float:
        """Get current memory usage in MB"""
        if self.device.type == 'cuda':
            torch.cuda.reset_peak_memory_stats(self.device)
            torch.cuda.synchronize(self.device)
            return torch.cuda.memory_allocated() / 1024 / 1024  # Convert to MB
        else:
            # For CPU, this is more complex and less accurate
            return 0.0
    
    def _create_inputs(self, ln_init, batch_size, seq_len):
        """
        Return a tensor of shape (B, N, N, D) drawn from N(beta, |gamma|^2)
        where beta and gamma are taken from `layer`.
        If D is None we infer it from layer.normalized_shape.
        """
        # --- sanity checks -------------------------------------------------------
        tail_shape = ln_init.normalized_shape  # tuple or int
        if isinstance(tail_shape, int):
            tail_shape = (tail_shape,)

        full_shape = (batch_size, seq_len, seq_len) + tail_shape           # (B, N, N, D)
        # --- broadcast beta and |gamma| over leading axes -----------------------
        beta = ln_init.bias.reshape((1,)*3 + tail_shape)
        gamma = ln_init.weight.abs().reshape((1,)*3 + tail_shape)

        # --- draw samples --------------------------------------------------------
        x = torch.normal(mean=beta.expand(full_shape),
                        std=gamma.expand(full_shape),
                        generator=None).to(device=self.device, dtype=self.precision)
        
        # Set requires_grad=True to enable gradient computation for backward pass testing
        x.requires_grad_(True)

        return x
    
    def _get_peak_memory_usage(self) -> float:
        """Get peak memory usage in MB"""
        if self.device.type == 'cuda':
            torch.cuda.synchronize(self.device)
            return torch.cuda.max_memory_allocated() / 1024 / 1024  # Convert to MB
        else:
            return 0.0

    def time_once(self, func, *args, **kwargs) -> float:
        start = torch.cuda.Event(enable_timing=True)
        end   = torch.cuda.Event(enable_timing=True)
        times = []
        
        # Warmup - create fresh inputs for each iteration
        for _ in range(self.num_warmup):
            # Create fresh arguments to avoid computational graph accumulation
            fresh_args = [arg.detach().clone() if isinstance(arg, torch.Tensor) and arg.requires_grad 
                         else arg for arg in args]
            _ = func(*fresh_args, **kwargs)
            if self.device.type == 'cuda':
                torch.cuda.synchronize()
            # Explicit cleanup
            torch.cuda.empty_cache()
        
        # Timing iterations with fresh inputs each time
        for _ in range(self.num_iterations):
            # Create fresh arguments to avoid computational graph accumulation  
            fresh_args = [arg.detach().clone() if isinstance(arg, torch.Tensor) and arg.requires_grad 
                         else arg for arg in args]
            
            start.record(stream=torch.cuda.current_stream(self.device))
            _ = func(*fresh_args, **kwargs)
            end.record(stream=torch.cuda.current_stream(self.device))
            end.synchronize()                           # waits here, not inside scope
            times.append(start.elapsed_time(end))       # milliseconds
            
            # Explicit cleanup after each timing iteration
            torch.cuda.empty_cache()
            
        return float(np.mean(times))

    def time_single_iteration(self, func, *args, **kwargs) -> float:
        """
        Simple single-iteration timing for quick comparison with full model timing.
        More realistic for comparing against actual training times.
        """
        # Single warmup
        fresh_args = [arg.detach().clone() if isinstance(arg, torch.Tensor) and arg.requires_grad 
                     else arg for arg in args]
        _ = func(*fresh_args, **kwargs)
        torch.cuda.synchronize()
        torch.cuda.empty_cache()
        
        # Single timed iteration
        fresh_args = [arg.detach().clone() if isinstance(arg, torch.Tensor) and arg.requires_grad 
                     else arg for arg in args]
        
        start = torch.cuda.Event(enable_timing=True)
        end = torch.cuda.Event(enable_timing=True)
        
        start.record(stream=torch.cuda.current_stream(self.device))
        _ = func(*fresh_args, **kwargs)
        end.record(stream=torch.cuda.current_stream(self.device))
        end.synchronize()
        
        return start.elapsed_time(end)  # milliseconds

    def run_once(self, func, *args, **kwargs) -> Tuple[float, float, torch.Tensor, torch.Tensor]:
        initial_memory = self._get_memory_usage()
        torch.cuda.reset_peak_memory_stats(self.device)
        
        # Create fresh arguments to avoid computational graph issues
        fresh_args = [arg.detach().clone() if isinstance(arg, torch.Tensor) and arg.requires_grad 
                     else arg for arg in args]
        
        result = func(*fresh_args, **kwargs)
        torch.cuda.synchronize(self.device)
        peak_alloc = torch.cuda.max_memory_allocated(self.device) / 1024 / 1024
        peak_rsvd  = torch.cuda.max_memory_reserved(self.device) / 1024 / 1024
        
        # Handle both single tensor return and (output, grad) tuple return
        if isinstance(result, tuple):
            output, grad = result
            return peak_alloc - initial_memory, peak_rsvd - initial_memory, output.float(), grad.float() if grad is not None else None
        else:
            # Legacy support for single output
            return peak_alloc - initial_memory, peak_rsvd - initial_memory, result.float(), None

    def test_backward_comparison(self):
        """
        Compare backward pass between different implementations.
        
        Note: Uses single-iteration timing by default for faster testing and more realistic 
        comparison with actual training times. Change time_single_iteration() to time_once() 
        for detailed multi-iteration benchmarks.
        """
        if self.device.type != 'cuda':
            self.skipTest("requires CUDA")
        
        print("\n" + "="*100)
        print("BACKWARD PERFORMANCE COMPARISON")
        print("="*100)
        print(f"{'Config':<45} {'Implementation':<15} {'Time (ms)':<12} {'Peak Mem (MB)':<15} {'Peak Reserved (MB)':<15}")
        print("-"*100)
        
        results = []
        
        with torch.amp.autocast(device_type=self.device.type, dtype=self.precision):
            for batch_size, seq_len, hidden_dim in self.train_test_configs:
                for direction in self.directions:
                    config_str = f"B={batch_size}, N={seq_len}, D={hidden_dim}, {direction}"
                    
                    # Test standard implementation
                    if direction == 'outgoing':
                        z = self._create_inputs(self.ln_init_outgoing, batch_size, seq_len)
                        standard_model = TriangleMultiplicationOutgoing(
                            c_z=hidden_dim, c_hidden=hidden_dim
                        )
                    else:
                        z = self._create_inputs(self.ln_init_incoming, batch_size, seq_len)
                        standard_model = TriangleMultiplicationIncoming(
                            c_z=hidden_dim, c_hidden=hidden_dim
                        )
                    
                    standard_model = standard_model.to(self.device)
                    
                    # Measure standard implementation
                    def standard_backward(z, mask):
                        # Create a fresh tensor with the same data but detached from any existing graph
                        z_input = z.detach().clone().requires_grad_(True)
                        standard_model.zero_grad()
                        output = standard_model(z_input, mask)
                        # Use retain_graph=True to allow multiple backward passes for timing
                        output.sum().backward(retain_graph=True)
                        return output, z_input.grad
                    
                    # Use single iteration timing for faster testing - change to time_once for detailed benchmarks
                    std_time = self.time_single_iteration(
                        standard_backward, z, mask=None
                    )
                    std_peak_mem, std_alloc, std_output, std_grad = self.run_once(
                        standard_backward, z, mask=None
                    )
                    
                    print(f"{config_str:<45} {'Standard':<15} {std_time:<12.3f} {std_peak_mem:<15.3f} {std_alloc:<15.3f}")
                    
                    # Test kernel implementation if available
                    kernel_available = False
                    fused_time = 0
                    fused_peak_mem = 0
                    fused_alloc = 0
                    fused_output = None

                    try:
                        def fused_backward(z, mask):
                            # Create a fresh tensor with the same data but detached from any existing graph
                            z_input = z.detach().clone().requires_grad_(True)
                            standard_model.zero_grad()
                            output = standard_model(z_input, mask, use_fused_kernel=True)
                            # Use retain_graph=True to allow multiple backward passes for timing
                            output.sum().backward(retain_graph=True)
                            return output, z_input.grad
                        
                        # Measure fused implementation  
                        fused_time = self.time_single_iteration(
                            fused_backward, z, mask=None
                        )
                        fused_peak_mem, fused_alloc, fused_output, fused_grad = self.run_once(
                            fused_backward, z, mask=None
                        )
                        
                        # Check forward pass consistency
                        self.assertTrue(
                            torch.allclose(std_output, fused_output, 
                                        atol=self.atol, rtol=self.rtol),
                            f"Forward outputs differ for config: batch={batch_size}, "
                            f"seq_len={seq_len}, hidden={hidden_dim}, direction={direction}\n"
                            f"Max diff: {(std_output - fused_output).abs().max().item()}"
                        )
                        
                        # Check backward pass (gradient) consistency
                        self.assertIsNotNone(std_grad, "Standard implementation should produce gradients")
                        self.assertIsNotNone(fused_grad, "Fused implementation should produce gradients") 
                        self.assertTrue(
                            torch.allclose(std_grad, fused_grad, 
                                        atol=self.atol, rtol=self.rtol),
                            f"Gradients differ for config: batch={batch_size}, "
                            f"seq_len={seq_len}, hidden={hidden_dim}, direction={direction}\n"
                            f"Max grad diff: {(std_grad - fused_grad).abs().max().item()}"
                        )
                        
                        kernel_available = True
                    except Exception as e:
                        warnings.warn(f"Kernel implementation not available: {e}")
                        kernel_available = False
                    
                    if kernel_available:
                        print(f"{'':<45} {'Fused':<15} {fused_time:<12.3f} {fused_peak_mem:<15.3f} {fused_alloc:<15.3f}")
                    
                        # Calculate improvements
                        time_improvement = (std_time - fused_time) / std_time * 100
                        mem_improvement = (std_peak_mem - fused_peak_mem) / std_peak_mem * 100
                        
                        print(f"{'':<45} {'Improvement':<15} {time_improvement:<12.1f}% {mem_improvement:<15.1f}%")
                        print("-"*100)
                    
                    results.append({
                        'config': config_str,
                        'std_time': std_time,
                        'std_mem': std_peak_mem,
                        'fused_time': fused_time,
                        'fused_mem': fused_peak_mem,
                        'time_improvement': time_improvement,
                        'mem_improvement': mem_improvement
                    })
        
        # Summary
        avg_time_improvement = np.mean([r['time_improvement'] for r in results])
        avg_mem_improvement = np.mean([r['mem_improvement'] for r in results])
        
        print(f"\nAVERAGE IMPROVEMENTS:")
        print(f"Time: {avg_time_improvement:.1f}% faster with fused implementation")
        print(f"Memory: {avg_mem_improvement:.1f}% less memory with fused implementation")
    
    def test_forward_comparison(self):
        """Compare memory usage between normal and inplace modes"""
        if self.device.type != 'cuda':
            self.skipTest("Memory measurement requires CUDA")
        
        print("\n" + "="*100)
        print("FORWARD PERFORMANCE COMPARISON")
        print("="*100)
        print(f"{'Config':<45} {'Mode':<15} {'Time (ms)':<12} {'Peak Mem (MB)':<15} {'Peak Reserved (MB)':<15}")
        print("-"*100)
        
        with torch.amp.autocast(device_type=self.device.type, dtype=self.precision):
            for batch_size, seq_len, hidden_dim in self.eval_test_configs:
                for direction in ['outgoing', 'incoming']:
                    config_str = f"B={batch_size}, N={seq_len}, D={hidden_dim}, {direction}"
                    
                    # Create inputs
                    mask = torch.ones(batch_size, seq_len, seq_len, 
                                    device=self.device, dtype=torch.float32)
                    
                    # Create model
                    if direction == 'outgoing':
                        z = self._create_inputs(self.ln_init_outgoing, batch_size, seq_len)
                        model = TriangleMultiplicationOutgoing(
                            c_z=hidden_dim, c_hidden=hidden_dim 
                        )
                    else:
                        z = self._create_inputs(self.ln_init_incoming, batch_size, seq_len)
                        model = TriangleMultiplicationIncoming(
                            c_z=hidden_dim, c_hidden=hidden_dim
                        )
                    
                    model = model.to(self.device).eval()
                    
                    # Measure normal mode
                    with torch.inference_mode():
                        normal_time = self.time_once(
                            model, z.clone(), mask=mask
                        )
                        normal_peak_mem, normal_alloc, normal_output, _ = self.run_once(
                            model, z.clone(), mask=mask
                        )
                    
                    print(f"{config_str:<45} {'Normal':<15} {normal_time:<12.3f} {normal_peak_mem:<15.3f} {normal_alloc:<15.3f}")
                    
                    # Measure fused mode
                    if seq_len < 2560:
                        with torch.inference_mode():
                            fused_time = self.time_once(
                                model, z.clone(), mask=mask, use_fused_kernel=True
                            )
                            fused_peak_mem, fused_alloc, fused_output, _ = self.run_once(
                                model, z.clone(), mask=mask, use_fused_kernel=True
                            )
                        
                        print(f"{'':<45} {'Fused':<15} {fused_time:<12.3f} {fused_peak_mem:<15.3f} {fused_alloc:<15.3f}")
                    
                    # Measure inplace mode
                    with torch.inference_mode():
                        def inplace_forward():
                            z_copy = z.clone()
                            return model(z_copy, mask=mask, inplace_safe=True, 
                                    _add_with_inplace=False, _inplace_chunk_size=128)
                        
                        inplace_time = self.time_once(inplace_forward)
                        inplace_peak_mem, inplace_alloc, inplace_output, _ = self.run_once(
                            inplace_forward
                        )
                    
                    print(f"{'':<45} {'Inplace':<15} {inplace_time:<12.3f} {inplace_peak_mem:<15.3f} {inplace_alloc:<15.3f}")
                    
                    # Measure fusion inplace mode
                    with torch.no_grad():
                        def inplace_forward():
                            z_copy = z.clone()
                            return model._fused_inference_forward(z_copy, mask=mask, 
                                    with_add=False, inplace_chunk_size=128)
                        
                        fused_inplace_time = self.time_once(inplace_forward)
                        fused_inplace_peak_mem, fused_inplace_alloc, fused_inplace_output, _ = self.run_once(
                            inplace_forward
                        )
                    
                    print(f"{'':<45} {'Fused Inplace':<15} {fused_inplace_time:<12.3f} {fused_inplace_peak_mem:<15.3f} {fused_inplace_alloc:<15.3f}")
                    
                    # check consistency
                    if seq_len < 2560:
                        self.assertTrue(
                            torch.allclose(normal_output, fused_output, 
                                        atol=self.atol, rtol=self.rtol),
                            f"Outputs differ for config: batch={batch_size}, "
                            f"seq_len={seq_len}, hidden={hidden_dim}, direction={direction}\n"
                            f"Max diff: {(normal_output - fused_output).abs().max().item()}"
                        )
                    
                    self.assertTrue(
                            torch.allclose(normal_output, inplace_output, 
                                        atol=self.atol, rtol=self.rtol),
                            f"Outputs differ for config: batch={batch_size}, "
                            f"seq_len={seq_len}, hidden={hidden_dim}, direction={direction}\n"
                            f"Max diff: {(normal_output - inplace_output).abs().max().item()}"
                        )
                    self.assertTrue(
                            torch.allclose(normal_output, fused_inplace_output, 
                                        atol=self.atol, rtol=self.rtol),
                            f"Outputs differ for config: batch={batch_size}, "
                            f"seq_len={seq_len}, hidden={hidden_dim}, direction={direction}\n"
                            f"Max diff: {(normal_output - fused_inplace_output).abs().max().item()}"
                        )
                    
                    # Calculate improvements
                    if seq_len < 2560:
                        time_diff = (fused_time - normal_time) / normal_time * 100
                        mem_improvement = (normal_peak_mem - fused_peak_mem) / normal_peak_mem * 100
                        print(f"{'':<45} {'Diff NvsF':<15} {time_diff:<12.1f}% {mem_improvement:<15.1f}% {(normal_output - fused_output).abs().max().item():.5f}")
                    
                    time_diff = (inplace_time - normal_time) / normal_time * 100
                    mem_improvement = (normal_peak_mem - inplace_peak_mem) / normal_peak_mem * 100
                    
                    print(f"{'':<45} {'Diff NvsI':<15} {time_diff:<12.1f}% {mem_improvement:<15.1f}% {(normal_output - inplace_output).abs().max().item():.5f}")
                    
                    time_diff = (fused_inplace_time - inplace_time) / inplace_time * 100
                    mem_improvement = (inplace_peak_mem - fused_inplace_peak_mem) / inplace_peak_mem * 100
                    print(f"{'':<45} {'Diff IvsFusedI':<15} {time_diff:<12.1f}% {mem_improvement:<15.1f}% {(inplace_output - fused_inplace_output).abs().max().item():.5f}")
                    
                    time_diff = (fused_inplace_time - normal_time) / normal_time * 100
                    mem_improvement = (normal_peak_mem - fused_inplace_peak_mem) / normal_peak_mem * 100
                    print(f"{'':<45} {'Diff NvsFusedI':<15} {time_diff:<12.1f}% {mem_improvement:<15.1f}% {(normal_output - fused_inplace_output).abs().max().item():.5f}")
                    
                    print("-"*100)


if __name__ == "__main__":
    unittest.main()
