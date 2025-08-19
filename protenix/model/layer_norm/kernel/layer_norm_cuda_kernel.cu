// Copyright 2024 ByteDance and/or its affiliates.
// Copyright 2020 The OneFlow Authors.
// Copyright 2021- HPC-AI Technology Inc.
//
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//
//      http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.
#include <cooperative_groups.h>
#include <cuda.h>
#include <cuda_runtime.h>
#include <torch/extension.h>
#include <iostream>

#include <THC/THCDeviceUtils.cuh>

#include "ATen/ATen.h"
#include "ATen/AccumulateType.h"
#include "ATen/cuda/CUDAContext.h"
#include "compat.h"
#include "type_shim.h"

#define CHECK_CUDA(x) TORCH_CHECK(x.is_cuda(), #x " must be a CUDA tensor")
#define CHECK_CONTIGUOUS(x) TORCH_CHECK(x.is_contiguous(), #x " must be contiguous")
#define CHECK_INPUT(x) \
    CHECK_CUDA(x);     \
    CHECK_CONTIGUOUS(x)

#define WarpNum 8
#define WarpSize 32
#define BlockSzie WarpNum*WarpSize

inline __device__ void WelfordOnline(float val, float* mean, float* m2, float* count) {
    *count += 1;
    float delta1 = val - *mean;
    *mean += delta1 / (*count);
    float delta2 = val - *mean;
    *m2 += delta1 * delta2;
}

inline __device__ void WelfordOnline(float b_mean, float b_m2, float b_count, float* mean,
                                     float* m2, float* count) {
    if (b_count == 0) {
        return;
    }
    float new_count = *count + b_count;
    float nb_n = b_count / new_count;
    float delta = b_mean - *mean;
    *mean += delta * nb_n;
    *m2 += b_m2 + delta * delta * (*count) * nb_n;
    *count = new_count;
}

__inline__ __device__ void WelfordWarpAllReduce(float thread_mean, float thread_m2,
                                                float thread_count, float* mean, float* m2,
                                                float* count, int syc_thread_num=32) {
    *mean = thread_mean;
    *m2 = thread_m2;
    *count = thread_count;
    for(int mask = syc_thread_num/2; mask >= 1; mask /= 2) {
        float b_mean = __shfl_xor_sync(0xffffffff, *mean, mask);
        float b_m2 = __shfl_xor_sync(0xffffffff, *m2, mask);
        float b_count = __shfl_xor_sync(0xffffffff, *count, mask);
        WelfordOnline(b_mean, b_m2, b_count, mean, m2, count);
    }
}

extern __shared__ float shared_data[];
template <typename T>
__global__ void LayerNormForward(T* input, T* output, T* gamma, T* beta, float* mean,
                                 float* invvar, int rows, int cols, double epsilon) {
    int warp_id = threadIdx.x / WarpSize;
    int lane_id = threadIdx.x % WarpSize;
    int row_offset = blockIdx.x * WarpNum + warp_id;

    T* shared_data_warp = (T*)shared_data + warp_id*cols;

    if (row_offset < rows) {
        T* row_input = input + (long long)(row_offset) * (long long)(cols); // Starting point for input data
        T* row_output = output + (long long)(row_offset) * (long long)(cols); // Starting point for output data

        float thread_mean = 0.f;
        float thread_m2 = 0.f;
        float thread_count = 0.f;

        float warp_mean;
        float warp_m2;
        float warp_count;
        // load data to shared memory
#pragma unroll
        for(int idx = lane_id; idx < cols; idx += WarpSize) {
            shared_data_warp[idx] = row_input[idx];
            WelfordOnline(static_cast<float>(shared_data_warp[idx]), &thread_mean, &thread_m2, &thread_count);
        }

        WelfordWarpAllReduce(thread_mean, thread_m2, thread_count, &warp_mean, &warp_m2,
                             &warp_count);

        float row_mean = warp_mean;
        float row_variance = max(warp_m2 / warp_count, 0.f);
        float row_inv_var = rsqrt(row_variance + epsilon);
        if (lane_id == 0) {
            mean[row_offset] = row_mean;
            invvar[row_offset] = row_inv_var;
        }
        int process_type = (gamma != NULL)*2 + (beta != NULL);
        if (process_type == 0) {
#pragma unroll
            for(int idx = lane_id; idx < cols; idx += WarpSize)
                row_output[idx] = static_cast<T>((static_cast<float>(shared_data_warp[idx]) - row_mean) * row_inv_var);
        } else if (process_type == 1) {
#pragma unroll
            for(int idx = lane_id; idx < cols; idx += WarpSize)
                row_output[idx] = static_cast<T>((static_cast<float>(shared_data_warp[idx]) - row_mean) * row_inv_var + beta[idx]);
        } else if(process_type == 2) {
#pragma unroll
            for(int idx = lane_id; idx < cols; idx += WarpSize)
                row_output[idx] = static_cast<T>((static_cast<float>(shared_data_warp[idx]) - row_mean) * row_inv_var * gamma[idx]);
        } else {
#pragma unroll
            for(int idx = lane_id; idx < cols; idx += WarpSize)
                row_output[idx] = static_cast<T>((static_cast<float>(shared_data_warp[idx]) - row_mean) * row_inv_var * gamma[idx] + beta[idx]);
        }
    }
}

int find_opt_threads(int M, int elements_per_thread) {
    const int candidates[] = {1, 2, 4, 8, 16, 32};
    int min_threads = (M + elements_per_thread - 1) / elements_per_thread; // 向上取整
    for (int k : candidates) {
        if (k >= min_threads) {
            return k;
        }
    }
    return 32;
}

template <typename T, typename VecType>
__global__ void LayerNormForwardV2(T* input, T* output, T* gamma, T* beta, 
                                   float* mean, float* invvar, long rows, 
                                   long cols, float epsilon) {
    constexpr int ELEMENTS_PER_THREAD = sizeof(VecType) / sizeof(T);
    
    const long tid = threadIdx.x;
    const long row = blockIdx.x * blockDim.y + threadIdx.y;
    T* row_input_ptr = input + row * cols;
    T* row_output_ptr = output + row * cols;
    if (row >= rows) return;
    
    const bool has_gamma = gamma != nullptr;
    const bool has_beta = beta != nullptr;

    const long ELEMENTS_PER_BLOCK = blockDim.x * ELEMENTS_PER_THREAD;
    long TOTAL_BLOCKS = cols / ELEMENTS_PER_BLOCK;

    {
        const long REMAINING_ELEMENTS = cols % ELEMENTS_PER_BLOCK;
        const long REMAINING_VECTORS = REMAINING_ELEMENTS / ELEMENTS_PER_THREAD;
        if(tid < REMAINING_VECTORS) ++TOTAL_BLOCKS;
    }

    float thread_mean = 0.f, thread_m2 = 0.f, thread_count = 0.f;
        for (int block = 0; block < TOTAL_BLOCKS; ++block) {
        const long base_idx = block * ELEMENTS_PER_BLOCK + tid * ELEMENTS_PER_THREAD;
        VecType vec = *reinterpret_cast<VecType*>(row_input_ptr + base_idx);
        T* vals = reinterpret_cast<T*>(&vec);
        #pragma unroll
        for (int i = 0; i < ELEMENTS_PER_THREAD; ++i) {
            WelfordOnline(static_cast<float>(vals[i]), &thread_mean, &thread_m2, &thread_count);
        }
    }

    float row_mean, warp_m2, warp_count;
    WelfordWarpAllReduce(thread_mean, thread_m2, thread_count, &row_mean, &warp_m2, &warp_count, blockDim.x);
    float row_inv_var = rsqrt(max(warp_m2 / warp_count, 0.f) + epsilon);

    if (tid == 0) {
        mean[row] = row_mean;
        invvar[row] = row_inv_var;
    }

    for (int block = 0; block < TOTAL_BLOCKS; ++block) {
        const long base_idx = block * ELEMENTS_PER_BLOCK + tid * ELEMENTS_PER_THREAD;
        VecType vec = *reinterpret_cast<VecType*>(row_input_ptr + base_idx);
        VecType vec_out;
        T* vals = reinterpret_cast<T*>(&vec);
        T* vals_out = reinterpret_cast<T*>(&vec_out);

        VecType vec_gamma;
        if (has_gamma) {
            vec_gamma = *reinterpret_cast<VecType*>(gamma + base_idx);
        } else {
            T default_gamma[ELEMENTS_PER_THREAD];
            for (int i = 0; i < ELEMENTS_PER_THREAD; ++i)
                default_gamma[i] = static_cast<T>(1.0f);
            vec_gamma = *reinterpret_cast<VecType*>(default_gamma);
        }

        VecType vec_beta;
        if (has_beta) {
            vec_beta = *reinterpret_cast<VecType*>(beta + base_idx);
        } else {
            T default_beta[ELEMENTS_PER_THREAD];
            for (int i = 0; i < ELEMENTS_PER_THREAD; ++i)
                default_beta[i] = static_cast<T>(0.0f);
            vec_beta = *reinterpret_cast<VecType*>(default_beta);
        }

        T* gamma_vals = reinterpret_cast<T*>(&vec_gamma);
        T* beta_vals = reinterpret_cast<T*>(&vec_beta);

        #pragma unroll
        for (int i = 0; i < ELEMENTS_PER_THREAD; ++i) {
            float normalized = (static_cast<float>(vals[i]) - row_mean) * row_inv_var;
            normalized = normalized * static_cast<float>(gamma_vals[i]) 
                       + static_cast<float>(beta_vals[i]);
            vals_out[i] = static_cast<T>(normalized);
        }

        *reinterpret_cast<VecType*>(row_output_ptr + base_idx) = vec_out;
    }
}

void cuda_layer_norm(at::Tensor* output, at::Tensor* mean, at::Tensor* invvar, at::Tensor* input,
                     int rows, int cols, at::IntArrayRef normalized_shape, at::Tensor* gamma,
                     at::Tensor* beta, double epsilon) {
    // 获取元素字节大小
    const auto dtype = output->dtype();
    int element_size;
    if (dtype == torch::kFloat32) {
        element_size = 4;
    } else if (dtype == torch::kFloat16 || dtype == torch::kBFloat16) {
        element_size = 2;
    } else {
        throw std::runtime_error("Unsupported data type");
    }

    // 计算总字节数并检查对齐
    const int total_bytes = cols * element_size;
    int vec_size = 0;
    if (total_bytes % 16 == 0) {
        vec_size = 16;
    } else if (total_bytes % 8 == 0) {
        vec_size = 8;
    } else if (total_bytes % 4 == 0) {
        vec_size = 4;
    } else {
        vec_size = 2;
    }

    // 计算每线程处理的元素数
    const int elements_per_thread = vec_size / element_size;
    const int threads_per_row = find_opt_threads(cols, elements_per_thread);

    // 配置kernel参数
    const int threads_per_block = 128;
    const int rows_per_block = threads_per_block / threads_per_row;
    const dim3 grid((rows + rows_per_block - 1) / rows_per_block);
    const dim3 block(threads_per_row, rows_per_block);

    // 类型分发
    if (dtype == torch::kFloat32) {
        if (vec_size == 16) {
            LayerNormForwardV2<float, float4><<<grid, block>>>(
                static_cast<float*>(input->data_ptr()),
                static_cast<float*>(output->data_ptr()),
                gamma ? static_cast<float*>(gamma->data_ptr()) : nullptr,
                beta ? static_cast<float*>(beta->data_ptr()) : nullptr,
                static_cast<float*>(mean->data_ptr()),
                static_cast<float*>(invvar->data_ptr()),
                long(rows), long(cols), float(epsilon)
            );
        } else if (vec_size == 8) {
            LayerNormForwardV2<float, float2><<<grid, block>>>(
                static_cast<float*>(input->data_ptr()),
                static_cast<float*>(output->data_ptr()),
                gamma ? static_cast<float*>(gamma->data_ptr()) : nullptr,
                beta ? static_cast<float*>(beta->data_ptr()) : nullptr,
                static_cast<float*>(mean->data_ptr()),
                static_cast<float*>(invvar->data_ptr()),
                long(rows), long(cols), float(epsilon)
            );
        } else if (vec_size == 4) {
            LayerNormForwardV2<float, float><<<grid, block>>>(
                static_cast<float*>(input->data_ptr()),
                static_cast<float*>(output->data_ptr()),
                gamma ? static_cast<float*>(gamma->data_ptr()) : nullptr,
                beta ? static_cast<float*>(beta->data_ptr()) : nullptr,
                static_cast<float*>(mean->data_ptr()),
                static_cast<float*>(invvar->data_ptr()),
                long(rows), long(cols), float(epsilon)
            );
        }
    } 
    else if (dtype == torch::kFloat16) {
                if (vec_size == 16) {  // 使用float4处理half类型（8个元素）
            LayerNormForwardV2<at::Half, float4><<<grid, block>>>(
                static_cast<at::Half*>(input->data_ptr()),
                static_cast<at::Half*>(output->data_ptr()),
                gamma ? static_cast<at::Half*>(gamma->data_ptr()) : nullptr,
                beta ? static_cast<at::Half*>(beta->data_ptr()) : nullptr,
                static_cast<float*>(mean->data_ptr()),
                static_cast<float*>(invvar->data_ptr()),
                long(rows), long(cols), float(epsilon)
            );
        } else if (vec_size == 8) {  // float2处理4个half元素
            LayerNormForwardV2<at::Half, float2><<<grid, block>>>(
                static_cast<at::Half*>(input->data_ptr()),
                static_cast<at::Half*>(output->data_ptr()),
                gamma ? static_cast<at::Half*>(gamma->data_ptr()) : nullptr,
                beta ? static_cast<at::Half*>(beta->data_ptr()) : nullptr,
                static_cast<float*>(mean->data_ptr()),
                static_cast<float*>(invvar->data_ptr()),
                long(rows), long(cols), float(epsilon)
            );
        } else if (vec_size == 4) {  // float处理2个half元素
            LayerNormForwardV2<at::Half, float><<<grid, block>>>(
                static_cast<at::Half*>(input->data_ptr()),
                static_cast<at::Half*>(output->data_ptr()),
                gamma ? static_cast<at::Half*>(gamma->data_ptr()) : nullptr,
                beta ? static_cast<at::Half*>(beta->data_ptr()) : nullptr,
                static_cast<float*>(mean->data_ptr()),
                static_cast<float*>(invvar->data_ptr()),
                long(rows), long(cols), float(epsilon)
            );
        } else if (vec_size == 2) {
            LayerNormForwardV2<at::Half, at::Half><<<grid, block>>>(
                static_cast<at::Half*>(input->data_ptr()),
                static_cast<at::Half*>(output->data_ptr()),
                gamma ? static_cast<at::Half*>(gamma->data_ptr()) : nullptr,
                beta ? static_cast<at::Half*>(beta->data_ptr()) : nullptr,
                static_cast<float*>(mean->data_ptr()),
                static_cast<float*>(invvar->data_ptr()),
                long(rows), long(cols), float(epsilon)
            );
        }
    }
    else if (dtype == torch::kBFloat16) {
        if (vec_size == 16) {
            LayerNormForwardV2<at::BFloat16, float4><<<grid, block>>>(
                static_cast<at::BFloat16*>(input->data_ptr()),
                static_cast<at::BFloat16*>(output->data_ptr()),
                gamma ? static_cast<at::BFloat16*>(gamma->data_ptr()) : nullptr,
                beta ? static_cast<at::BFloat16*>(beta->data_ptr()) : nullptr,
                static_cast<float*>(mean->data_ptr()),
                static_cast<float*>(invvar->data_ptr()),
                long(rows), long(cols), float(epsilon)
            );
        } else if (vec_size == 8) {
            LayerNormForwardV2<at::BFloat16, float2><<<grid, block>>>(
                static_cast<at::BFloat16*>(input->data_ptr()),
                static_cast<at::BFloat16*>(output->data_ptr()),
                gamma ? static_cast<at::BFloat16*>(gamma->data_ptr()) : nullptr,
                beta ? static_cast<at::BFloat16*>(beta->data_ptr()) : nullptr,
                static_cast<float*>(mean->data_ptr()),
                static_cast<float*>(invvar->data_ptr()),
                long(rows), long(cols), float(epsilon)
            );
        } else if (vec_size == 4) {
            LayerNormForwardV2<at::BFloat16, float><<<grid, block>>>(
                static_cast<at::BFloat16*>(input->data_ptr()),
                static_cast<at::BFloat16*>(output->data_ptr()),
                gamma ? static_cast<at::BFloat16*>(gamma->data_ptr()) : nullptr,
                beta ? static_cast<at::BFloat16*>(beta->data_ptr()) : nullptr,
                static_cast<float*>(mean->data_ptr()),
                static_cast<float*>(invvar->data_ptr()),
                long(rows), long(cols), float(epsilon)
            );
        } else if (vec_size == 2) {
            LayerNormForwardV2<at::BFloat16, at::BFloat16><<<grid, block>>>(
                static_cast<at::BFloat16*>(input->data_ptr()),
                static_cast<at::BFloat16*>(output->data_ptr()),
                gamma ? static_cast<at::BFloat16*>(gamma->data_ptr()) : nullptr,
                beta ? static_cast<at::BFloat16*>(beta->data_ptr()) : nullptr,
                static_cast<float*>(mean->data_ptr()),
                static_cast<float*>(invvar->data_ptr()),
                long(rows), long(cols), float(epsilon)
            );
        }
    } 
}

template <typename T>
struct SharedMemory;

template <>
struct SharedMemory<float> {
    __device__ float* getPointer() {
        extern __shared__ float s_float[];
        return s_float;
    }
};

template<typename T>
__inline__ __device__ T WarpReduce(T val) {
  for (int mask = 16; mask > 0; mask /= 2) { val += __shfl_xor_sync(0xffffffff, val, mask); }
  return val;
}

constexpr int tile_size = 32;
constexpr int num_per_block = 4;
constexpr int block_dim_x = 32;
constexpr int block_dim_y = 32 / num_per_block;

template <typename T, typename V>
__global__ void LayerNormParamGradStep1(int rows, int cols, const V* __restrict__ dy,
                                        const T* __restrict__ x, const float* __restrict__ mean,
                                        const float* __restrict__ inv_var,
                                        float* __restrict__ tmp_gamma_diff, float* __restrict__ tmp_beta_diff) {
  __shared__ float dgamma[32][33];
  __shared__ float dbeta[32][33];
  float dgamma_sum[num_per_block];
  float dbeta_sum[num_per_block];
#pragma unroll
  for (int index = 0; index < num_per_block; ++index) {
    dgamma_sum[index] = 0;
    dbeta_sum[index] = 0;
  }
  const int col_id = blockIdx.x * blockDim.x + threadIdx.x;
  if (col_id < cols) {
    for (int i = blockIdx.y * tile_size + threadIdx.y; i < rows; i += tile_size * gridDim.y) {
#pragma unroll
      for (int index = 0; index < num_per_block; ++index) {
        int row_id = i + index * blockDim.y;
        if (row_id < rows) {
          int offset = row_id * cols + col_id;
          const float dy_val = static_cast<float>(dy[offset]);
          const float x_val = static_cast<float>(x[offset]);
          const float mean_val = mean[row_id];
          const float inv_var_val = inv_var[row_id];
          dgamma_sum[index] += dy_val * (x_val - mean_val) * inv_var_val;
          dbeta_sum[index] += dy_val;
        }
      }
    }
  }
#pragma unroll
  for (int index = 0; index < num_per_block; ++index) {
    dgamma[index * blockDim.y + threadIdx.y][threadIdx.x] = dgamma_sum[index];
    dbeta[index * blockDim.y + threadIdx.y][threadIdx.x] = dbeta_sum[index];
  }
  __syncthreads();
#pragma unroll
  for (int index = 0; index < num_per_block; ++index) {
    const int col_id = blockIdx.x * blockDim.x + threadIdx.y + index * blockDim.y;
    if (col_id < cols) {
      float gamma_sum = dgamma[threadIdx.x][threadIdx.y + index * blockDim.y];
      float beta_sum = dbeta[threadIdx.x][threadIdx.y + index * blockDim.y];
      float global_dgamma = WarpReduce<float>(gamma_sum);
      float global_dbeta = WarpReduce<float>(beta_sum);
      if (threadIdx.x == 0) {
        const int offset = blockIdx.y * cols + col_id;
        tmp_gamma_diff[offset] = global_dgamma;
        tmp_beta_diff[offset] = global_dbeta;
      }
    }
  }
}

template <typename T, typename V>
__global__ void LayerNormGammaGradStep1(int rows, int cols, const V* __restrict__ dy,
                                        const T* __restrict__ x, const float* __restrict__ mean,
                                        const float* __restrict__ inv_var, float* __restrict__ tmp_gamma_diff) {
  __shared__ float dgamma[32][33];
  float dgamma_sum[num_per_block];
#pragma unroll
  for (int index = 0; index < num_per_block; ++index) {
    dgamma_sum[index] = 0;
  }
  const int col_id = blockIdx.x * blockDim.x + threadIdx.x;
  if (col_id < cols) {
    for (int i = blockIdx.y * tile_size + threadIdx.y; i < rows; i += tile_size * gridDim.y) {
#pragma unroll
      for (int index = 0; index < num_per_block; ++index) {
        int row_id = i + index * blockDim.y;
        if (row_id < rows) {
          int offset = row_id * cols + col_id;
          const float dy_val = static_cast<float>(dy[offset]);
          const float x_val = static_cast<float>(x[offset]);
          const float mean_val = mean[row_id];
          const float inv_var_val = inv_var[row_id];
          dgamma_sum[index] += dy_val * (x_val - mean_val) * inv_var_val;
        }
      }
    }
  }
#pragma unroll
  for (int index = 0; index < num_per_block; ++index) {
    dgamma[index * blockDim.y + threadIdx.y][threadIdx.x] = dgamma_sum[index];
  }
  __syncthreads();
#pragma unroll
  for (int index = 0; index < num_per_block; ++index) {
    const int col_id = blockIdx.x * blockDim.x + threadIdx.y + index * blockDim.y;
    if (col_id < cols) {
      float gamma_sum = dgamma[threadIdx.x][threadIdx.y + index * blockDim.y];
      float global_dgamma = WarpReduce<float>(gamma_sum);
      if (threadIdx.x == 0) {
        const int offset = blockIdx.y * cols + col_id;
        tmp_gamma_diff[offset] = global_dgamma;
      }
    }
  }
}

template <typename T, typename V>
__global__ void LayerNormBetaGradStep1(int rows, int cols, const V* __restrict__ dy,
                                        const T* __restrict__ x, const float* __restrict__ mean,
                                        const float* __restrict__ inv_var, float* __restrict__ tmp_beta_diff) {
  __shared__ float dbeta[32][33];
  float dbeta_sum[num_per_block];
#pragma unroll
  for (int index = 0; index < num_per_block; ++index) {
    dbeta_sum[index] = 0;
  }
  const int col_id = blockIdx.x * blockDim.x + threadIdx.x;
  if (col_id < cols) {
    for (int i = blockIdx.y * tile_size + threadIdx.y; i < rows; i += tile_size * gridDim.y) {
#pragma unroll
      for (int index = 0; index < num_per_block; ++index) {
        int row_id = i + index * blockDim.y;
        if (row_id < rows) {
          int offset = row_id * cols + col_id;
          const float dy_val = static_cast<float>(dy[offset]);
          dbeta_sum[index] += dy_val;
        }
      }
    }
  }
#pragma unroll
  for (int index = 0; index < num_per_block; ++index) {
    dbeta[index * blockDim.y + threadIdx.y][threadIdx.x] = dbeta_sum[index];
  }
  __syncthreads();
#pragma unroll
  for (int index = 0; index < num_per_block; ++index) {
    const int col_id = blockIdx.x * blockDim.x + threadIdx.y + index * blockDim.y;
    if (col_id < cols) {
      float beta_sum = dbeta[threadIdx.x][threadIdx.y + index * blockDim.y];
      float global_dbeta = WarpReduce<float>(beta_sum);
      if (threadIdx.x == 0) {
        const int offset = blockIdx.y * cols + col_id;
        tmp_beta_diff[offset] = global_dbeta;
      }
    }
  }
}

template <typename V>
__global__ void LayerNormParamGradStep2(const float* part_grad_gamma, const float* part_grad_beta,
                                        const int part_size, const int row, const int col,
                                        V* grad_gamma, V* grad_beta) {
    // sum partial gradients for gamma and beta
    SharedMemory<float> shared;
    float* buf = shared.getPointer();
    int i2 = blockIdx.x * blockDim.x + threadIdx.x;
    if (i2 < col) {
        // each warp does sequential reductions until reduced part_size is num_warps
        // int num_warp_reductions = part_size / blockDim.y;
        float sum_gamma = float(0);
        float sum_beta = float(0);
        const float* part_grad_gamma_ptr = part_grad_gamma + i2;
        const float* part_grad_beta_ptr = part_grad_beta + i2;
        for (int row_idx = threadIdx.y; row_idx < part_size; row_idx += blockDim.y) {
            sum_gamma += part_grad_gamma_ptr[row_idx * col];
            sum_beta += part_grad_beta_ptr[row_idx * col];
        }
        // inter-warp reductions
        const int nbsize3 = blockDim.x * blockDim.y / 2;
        for (int offset = blockDim.y / 2; offset >= 1; offset /= 2) {
            // top half write to shared memory
            if (threadIdx.y >= offset && threadIdx.y < 2 * offset) {
                const int write_idx = (threadIdx.y - offset) * blockDim.x + threadIdx.x;
                buf[write_idx] = sum_gamma;
                buf[write_idx + nbsize3] = sum_beta;
            }
            __syncthreads();
            // bottom half sums
            if (threadIdx.y < offset) {
                const int read_idx = threadIdx.y * blockDim.x + threadIdx.x;
                sum_gamma += buf[read_idx];
                sum_beta += buf[read_idx + nbsize3];
            }
            __syncthreads();
        }
        // write out fully summed gradients
        if (threadIdx.y == 0) {
            grad_gamma[i2] = sum_gamma;
            grad_beta[i2] = sum_beta;
        }
    }
}

template <typename V>
__global__ void LayerNormGammaGradStep2(const float* part_grad_gamma, const int part_size, const int row, const int col, V* grad_gamma) {
    // sum partial gradients for gamma and beta
    SharedMemory<float> shared;
    float* buf = shared.getPointer();
    int i2 = blockIdx.x * blockDim.x + threadIdx.x;
    if (i2 < col) {
        // each warp does sequential reductions until reduced part_size is num_warps
        // int num_warp_reductions = part_size / blockDim.y;
        float sum_gamma = float(0);
        const float* part_grad_gamma_ptr = part_grad_gamma + i2;
        for (int row_idx = threadIdx.y; row_idx < part_size; row_idx += blockDim.y) {
            sum_gamma += part_grad_gamma_ptr[row_idx * col];
        }
        // inter-warp reductions
        for (int offset = blockDim.y / 2; offset >= 1; offset /= 2) {
            // top half write to shared memory
            if (threadIdx.y >= offset && threadIdx.y < 2 * offset) {
                const int write_idx = (threadIdx.y - offset) * blockDim.x + threadIdx.x;
                buf[write_idx] = sum_gamma;
            }
            __syncthreads();
            // bottom half sums
            if (threadIdx.y < offset) {
                const int read_idx = threadIdx.y * blockDim.x + threadIdx.x;
                sum_gamma += buf[read_idx];
            }
            __syncthreads();
        }
        // write out fully summed gradients
        if (threadIdx.y == 0) {
            grad_gamma[i2] = sum_gamma;
        }
    }
}

template <typename V>
__global__ void LayerNormBetaGradStep2(const float* part_grad_beta, const int part_size, const int row, const int col, V* grad_beta) {
    // sum partial gradients for gamma and beta
    SharedMemory<float> shared;
    float* buf = shared.getPointer();
    int i2 = blockIdx.x * blockDim.x + threadIdx.x;
    if (i2 < col) {
        // each warp does sequential reductions until reduced part_size is num_warps
        // int num_warp_reductions = part_size / blockDim.y;
        float sum_beta = float(0);
        const float* part_grad_beta_ptr = part_grad_beta + i2;
        for (int row_idx = threadIdx.y; row_idx < part_size; row_idx += blockDim.y) {
            sum_beta += part_grad_beta_ptr[row_idx * col];
        }
        // inter-warp reductions
        for (int offset = blockDim.y / 2; offset >= 1; offset /= 2) {
            // top half write to shared memory
            if (threadIdx.y >= offset && threadIdx.y < 2 * offset) {
                const int write_idx = (threadIdx.y - offset) * blockDim.x + threadIdx.x;
                buf[write_idx] = sum_beta;
            }
            __syncthreads();
            // bottom half sums
            if (threadIdx.y < offset) {
                const int read_idx = threadIdx.y * blockDim.x + threadIdx.x;
                sum_beta += buf[read_idx];
            }
            __syncthreads();
        }
        // write out fully summed gradients
        if (threadIdx.y == 0) {
            grad_beta[i2] = sum_beta;
        }
    }
}

template <typename T, typename V>
__global__ void LayerNormInputGrad(const V* __restrict__ dout, const T* __restrict__ input,
                                   const int rows, const int cols, const float* __restrict__ mean,
                                   const float* __restrict__ invvar, float epsilon, const V* gamma,
                                   T* grad_input) {
    int WarpPerBlock = blockDim.x / WarpSize;
    int thread_idx = threadIdx.x;
    int warp_idx = thread_idx / WarpSize;
    int lane_idx = thread_idx % WarpSize;

    float* shared_dout = shared_data + warp_idx*cols;
    float* shared_input = shared_data + WarpPerBlock*cols + warp_idx*cols;
    float* shared_gamma = shared_data + 2*WarpPerBlock*cols;
    int row_stride = gridDim.x*WarpPerBlock;
    for(int row = blockIdx.x*WarpPerBlock+warp_idx; row < rows; row += row_stride) {
        float mean_r = mean[row];
        float invvar_r = invvar[row];
        // load dout, input and gamma
        long long data_offset = (long long)(row) * cols;
        const V* dout_r = dout + data_offset;
        const T* input_r = input + data_offset;
        T* grad_input_r = grad_input + data_offset;
#pragma unroll
        for(int col = lane_idx; col < cols; col += WarpSize) {
            shared_dout[col] = float(dout_r[col]);
            shared_input[col] = float(input_r[col]);
        }
        if(warp_idx == 0) {
#pragma unroll
            for(int col = lane_idx; col < cols; col += WarpSize) {
                shared_gamma[col] = gamma != NULL ? float(gamma[col]) : 1.0f;
            }
        }
        __syncthreads();

        float gamma_dout = 0.0;
        float gamma_dout_input_mean = 0.0;
        // reduction, gamma*dout and gamma*dout*(input-mean)
#pragma unroll
        for(int col = lane_idx; col < cols; col += WarpSize) {
            float temp = shared_gamma[col] * shared_dout[col];
            gamma_dout += temp;
            gamma_dout_input_mean += temp * (shared_input[col] - mean_r);
        }
        float global_gamma_dout = WarpReduce<float>(gamma_dout);
        float global_gamma_dout_input_mean = WarpReduce<float>(gamma_dout_input_mean);

        float part3_temp_value = global_gamma_dout_input_mean * invvar_r * invvar_r * invvar_r / cols;
        float part2 = global_gamma_dout * invvar_r / cols;
#pragma unroll
        for(int col = lane_idx; col < cols; col += WarpSize) {
            float part1 = shared_gamma[col] * shared_dout[col] * invvar_r;
            float part3 = (shared_input[col] - mean_r) * part3_temp_value;
            grad_input_r[col] = part1 - part2 - part3;
        }
    }
}

__inline__ __device__ void warp_sum_reduce(float &val, int syc_thread_num) {
    for(int mask = syc_thread_num/2; mask >= 1; mask /= 2) {
        val += __shfl_xor_sync(0xffffffff, val, mask);
    }
}

template <typename T, typename VecType>
__global__ void LayerNormInputGradV2(T* __restrict__ grad_output,
                                     T* __restrict__ input,
                                     int rows, int cols,
                                     float* __restrict__ mean,
                                     float* __restrict__ invvar,
                                     float epsilon, T* gamma,
                                     T* grad_input) {
    constexpr int ELEMENTS_PER_THREAD = sizeof(VecType) / sizeof(T);
    const int tid = threadIdx.x;
    const int row = blockIdx.x * blockDim.y + threadIdx.y;
    if (row >= rows) return;
    const int ELEMENTS_PER_BLOCK = blockDim.x * ELEMENTS_PER_THREAD;
    int TOTAL_BLOCKS = cols / ELEMENTS_PER_BLOCK;
    const int REMAINING_ELEMENTS = cols % ELEMENTS_PER_BLOCK;
    const int REMAINING_VECTORS = REMAINING_ELEMENTS / ELEMENTS_PER_THREAD;
    if (tid < REMAINING_VECTORS) ++TOTAL_BLOCKS;

    const float mean_val = mean[row];
    const float invvar_val = invvar[row];

    T* grad_output_row = grad_output + row * cols;
    T* input_row = input + row * cols;

    float gamma_mul_grad_output = 0.0;
    float gamma_mul_grad_output_input_mean = 0.0;
    const bool has_gamma = gamma != nullptr;
    for (int block = 0; block < TOTAL_BLOCKS; ++block) {
        const int base_idx = block * ELEMENTS_PER_BLOCK + tid * ELEMENTS_PER_THREAD;

        VecType gamma_vec;
        if (has_gamma) {
            gamma_vec = *reinterpret_cast<VecType*>(gamma + base_idx);
        } else {
            T default_gamma[ELEMENTS_PER_THREAD];
            for (int i = 0; i < ELEMENTS_PER_THREAD; ++i)
                default_gamma[i] = static_cast<T>(1.0f);
            gamma_vec = *reinterpret_cast<VecType*>(default_gamma);
        }
        VecType grad_output_vec = *reinterpret_cast<VecType*>(grad_output_row + base_idx);
        VecType input_vec = *reinterpret_cast<VecType*>(input_row + base_idx);

        T* gamma_vals = reinterpret_cast<T*>(&gamma_vec);
        T* grad_output_vals = reinterpret_cast<T*>(&grad_output_vec);
        
        #pragma unroll
        for (int i = 0; i < ELEMENTS_PER_THREAD; ++i) {
            gamma_mul_grad_output +=  gamma_vals[i] * grad_output_vals[i];
}

        T* input_vals = reinterpret_cast<T*>(&input_vec);
        #pragma unroll
        for (int i = 0; i < ELEMENTS_PER_THREAD; ++i) {
            gamma_mul_grad_output_input_mean += gamma_vals[i] * grad_output_vals[i] * (input_vals[i] - mean_val);
        }
    }
    warp_sum_reduce(gamma_mul_grad_output, blockDim.x);
    warp_sum_reduce(gamma_mul_grad_output_input_mean, blockDim.x);
    
    // 阶段2：计算公共系数
    const float k1 = gamma_mul_grad_output * invvar_val / cols;
    const float k2 = gamma_mul_grad_output_input_mean * invvar_val * invvar_val * invvar_val / cols;

    T* grad_input_row = grad_input + row * cols;
    // 阶段3：向量化写回梯度
    for (int block = 0; block < TOTAL_BLOCKS; ++block) {
        const int base_idx = block * ELEMENTS_PER_BLOCK + tid * ELEMENTS_PER_THREAD;

        // 重新加载必要数据
        VecType grad_vec = *reinterpret_cast<const VecType*>(grad_output_row + base_idx);
        VecType input_vec = *reinterpret_cast<const VecType*>(input_row + base_idx);
        VecType gamma_vec;
        if (has_gamma) {
            gamma_vec = *reinterpret_cast<VecType*>(gamma + base_idx);
        } else {
            T default_gamma[ELEMENTS_PER_THREAD];
            for (int i = 0; i < ELEMENTS_PER_THREAD; ++i)
                default_gamma[i] = static_cast<T>(1.0f);
            gamma_vec = *reinterpret_cast<VecType*>(default_gamma);
        }
        T* grad_vals = reinterpret_cast<T*>(&grad_vec);
        T* input_vals = reinterpret_cast<T*>(&input_vec);
        T* gamma_vals = reinterpret_cast<T*>(&gamma_vec);

        // 计算梯度
        VecType grad_input_vec;
        T* grad_input_vals = reinterpret_cast<T*>(&grad_input_vec);

        #pragma unroll
        for (int i = 0; i < ELEMENTS_PER_THREAD; ++i) {
            const float grad_val = static_cast<float>(grad_vals[i]);
            const float input_val = static_cast<float>(input_vals[i]);
            const float gamma_val = static_cast<float>(gamma_vals[i]);

            float grad = gamma_val * grad_val * invvar_val;
            grad -= k1;
            grad -= (input_val - mean_val) * k2;
            grad_input_vals[i] = static_cast<T>(grad);
        }

        // 向量化存储
        *reinterpret_cast<VecType*>(grad_input_row + base_idx) = grad_input_vec;
    }
}


template <typename T, typename V>
int GetGirdDimY(const int64_t num_instances, const int64_t norm_size) {
    const int grid_dim_x = (norm_size + tile_size - 1) / tile_size;
    const int max_grid_dim_y = (num_instances + tile_size - 1) / tile_size;
    const int block_size = block_dim_x * block_dim_y;
    int max_active_blocks = 0;
    cudaOccupancyMaxActiveBlocksPerMultiprocessor(
        &max_active_blocks, LayerNormParamGradStep1<T, V>, block_size, 0);
    int waves = 1;
    int dev;
    cudaGetDevice(&dev);
    int sm_count;
    cudaDeviceGetAttribute(&sm_count, cudaDevAttrMultiProcessorCount, dev);
    int num_blocks = max_active_blocks * sm_count * waves;
    int grid_dim_y = std::min(max_grid_dim_y, static_cast<int>(num_blocks / grid_dim_x));
    return std::max(grid_dim_y, 1);
}

template <typename T, typename V>
void HostLayerNormGradient(const V* dout, const float* mean, const float* invvar, at::Tensor* input, int row,
                           int col, const V* gamma, const V* beta, double epsilon, T* grad_input,
                           V* grad_gamma, V* grad_beta) {
    auto stream = at::cuda::getCurrentCUDAStream().stream();

    if (gamma != NULL && beta != NULL) {
        // compute grad_gamma(j) and grad_beta(j)
        const int part_size = GetGirdDimY<T, V>(row, col);
        const int grid_dim_x = (col + tile_size - 1) / tile_size;
        const int grid_dim_y = part_size;

        at::Tensor part_grad_gamma = at::empty({part_size, col}, input->options().dtype(at::ScalarType::Float));
        at::Tensor part_grad_beta = at::empty_like(part_grad_gamma);
        LayerNormParamGradStep1<T, V><<<dim3(grid_dim_x, grid_dim_y), dim3(32, 32 / num_per_block)>>>(
            row, col, dout, input->DATA_PTR<T>(), mean, invvar, part_grad_gamma.DATA_PTR<float>(), part_grad_beta.DATA_PTR<float>()
        );

        const dim3 threads3(32, 8, 1);
        const dim3 blocks3((col + 32 - 1) / 32, 1, 1);
        const int nshared3 = threads3.x * threads3.y * sizeof(float);
        LayerNormParamGradStep2<<<blocks3, threads3, nshared3, stream>>>(
            part_grad_gamma.DATA_PTR<float>(), part_grad_beta.DATA_PTR<float>(), part_size, row, col,
            grad_gamma, grad_beta);
    } else if (gamma != NULL && beta == NULL) {
        // compute grad_gamma(j) and grad_beta(j)
        const int part_size = GetGirdDimY<T, V>(row, col);
        const int grid_dim_x = (col + tile_size - 1) / tile_size;
        const int grid_dim_y = part_size;

        at::Tensor part_grad_gamma = at::empty({part_size, col}, input->options().dtype(at::ScalarType::Float));
        LayerNormGammaGradStep1<T, V><<<dim3(grid_dim_x, grid_dim_y), dim3(32, 32 / num_per_block)>>>(
            row, col, dout, input->DATA_PTR<T>(), mean, invvar, part_grad_gamma.DATA_PTR<float>());

        const dim3 threads3(32, 8, 1);
        const dim3 blocks3((col + 32 - 1) / 32, 1, 1);
        const int nshared3 = threads3.x * threads3.y * sizeof(float);
        LayerNormGammaGradStep2<<<blocks3, threads3, nshared3, stream>>>(
            part_grad_gamma.DATA_PTR<float>(), part_size, row, col, grad_gamma);
    } else if (gamma == NULL && beta!= NULL) {
        // compute grad_gamma(j) and grad_beta(j)
        const int part_size = GetGirdDimY<T, V>(row, col);
        const int grid_dim_x = (col + tile_size - 1) / tile_size;
        const int grid_dim_y = part_size;

        at::Tensor part_grad_beta = at::empty({part_size, col}, input->options().dtype(at::ScalarType::Float));
        LayerNormBetaGradStep1<T, V><<<dim3(grid_dim_x, grid_dim_y), dim3(32, 32 / num_per_block)>>>(
            row, col, dout, input->DATA_PTR<T>(), mean, invvar, part_grad_beta.DATA_PTR<float>()
        );

        const dim3 threads3(32, 8, 1);
        const dim3 blocks3((col + 32 - 1) / 32, 1, 1);
        const int nshared3 = threads3.x * threads3.y * sizeof(float);
        LayerNormBetaGradStep2<<<blocks3, threads3, nshared3, stream>>>(
            part_grad_beta.DATA_PTR<float>(), part_size, row, col, grad_beta);
    }

    // 获取元素字节大小
    const auto dtype = input->dtype();
    int element_size;
    if (dtype == torch::kFloat32) {
        element_size = 4;
    } else if (dtype == torch::kFloat16 || dtype == torch::kBFloat16) {
        element_size = 2;
    } else {
        throw std::runtime_error("Unsupported data type");
    }

    // 计算总字节数并检查对齐
    const int total_bytes = col * element_size;
    int vec_size = 0;
    if (total_bytes % 16 == 0) {
        vec_size = 16;
    } else if (total_bytes % 8 == 0) {
        vec_size = 8;
    } else if (total_bytes % 4 == 0) {
        vec_size = 4;
    } else {
        vec_size = 2;
    }

    // 计算每线程处理的元素数
    const int elements_per_thread = vec_size / element_size;
    const int threads_per_row = find_opt_threads(col, elements_per_thread);

    // 配置kernel参数
    const int threads_per_block = 128;
    const int rows_per_block = threads_per_block / threads_per_row;
    const dim3 grid((row + rows_per_block - 1) / rows_per_block);
    const dim3 block(threads_per_row, rows_per_block);
    if (dtype == torch::kFloat32) {
    if (vec_size == 16)
        LayerNormInputGradV2<float, float4><<<grid, block>>>((float*)dout, input->DATA_PTR<float>(), row, col, (float*)mean, (float*)invvar, float(epsilon), (float*)gamma, (float*)grad_input);
    else if(vec_size == 8)
        LayerNormInputGradV2<float, float2><<<grid, block>>>((float*)dout, input->DATA_PTR<float>(), row, col, (float*)mean, (float*)invvar, float(epsilon), (float*)gamma, (float*)grad_input);
    else if(vec_size == 4)
        LayerNormInputGradV2<float, float><<<grid, block>>>((float*)dout, input->DATA_PTR<float>(), row, col, (float*)mean, (float*)invvar, float(epsilon), (float*)gamma, (float*)grad_input);
    } else if (dtype == torch::kFloat16) {
    if (vec_size == 16)
        LayerNormInputGradV2<at::Half, float4><<<grid, block>>>((at::Half*)dout, input->DATA_PTR<at::Half>(), row, col, (float*)mean, (float*)invvar, float(epsilon), (at::Half*)gamma, (at::Half*)grad_input);
    else if(vec_size == 8)
        LayerNormInputGradV2<at::Half, float2><<<grid, block>>>((at::Half*)dout, input->DATA_PTR<at::Half>(), row, col, (float*)mean, (float*)invvar, float(epsilon), (at::Half*)gamma, (at::Half*)grad_input);
    else if(vec_size == 4)
        LayerNormInputGradV2<at::Half, float><<<grid, block>>>((at::Half*)dout, input->DATA_PTR<at::Half>(), row, col, (float*)mean, (float*)invvar, float(epsilon), (at::Half*)gamma, (at::Half*)grad_input);
    else if(vec_size == 2)
        LayerNormInputGradV2<at::Half, at::Half><<<grid, block>>>((at::Half*)dout, input->DATA_PTR<at::Half>(), row, col, (float*)mean, (float*)invvar, float(epsilon), (at::Half*)gamma, (at::Half*)grad_input);
    } else {
    if (vec_size == 16)
        LayerNormInputGradV2<at::BFloat16, float4><<<grid, block>>>((at::BFloat16*)dout, input->DATA_PTR<at::BFloat16>(), row, col, (float*)mean, (float*)invvar, float(epsilon), (at::BFloat16*)gamma, (at::BFloat16*)grad_input);
    else if(vec_size == 8)
        LayerNormInputGradV2<at::BFloat16, float2><<<grid, block>>>((at::BFloat16*)dout, input->DATA_PTR<at::BFloat16>(), row, col, (float*)mean, (float*)invvar, float(epsilon), (at::BFloat16*)gamma, (at::BFloat16*)grad_input);
    else if(vec_size == 4)
        LayerNormInputGradV2<at::BFloat16, float><<<grid, block>>>((at::BFloat16*)dout, input->DATA_PTR<at::BFloat16>(), row, col, (float*)mean, (float*)invvar, float(epsilon), (at::BFloat16*)gamma, (at::BFloat16*)grad_input);
    else if(vec_size == 2)
        LayerNormInputGradV2<at::BFloat16, at::BFloat16><<<grid, block>>>((at::BFloat16*)dout, input->DATA_PTR<at::BFloat16>(), row, col, (float*)mean, (float*)invvar, float(epsilon), (at::BFloat16*)gamma, (at::BFloat16*)grad_input);
    }
}

void cuda_layer_norm_gradient(at::Tensor* dout, at::Tensor* mean, at::Tensor* invvar,
                              at::Tensor* input, int row, int col, at::IntArrayRef normalized_shape,
                              at::Tensor* gamma, at::Tensor* beta, double epsilon,
                              at::Tensor* grad_input, at::Tensor* grad_gamma,
                              at::Tensor* grad_beta) {
    using namespace at;
    DISPATCH_FLOAT_HALF_AND_BFLOAT_INOUT_TYPES(
        input->scalar_type(), dout->scalar_type(), "cuda_layer_norm_gradient_kernel",
        HostLayerNormGradient(dout->DATA_PTR<scalar_t_out>(), mean->DATA_PTR<float>(),
                              invvar->DATA_PTR<float>(), input, row, col,
                              gamma != NULL ? gamma->DATA_PTR<scalar_t_out>() : NULL,
                              beta != NULL ? beta->DATA_PTR<scalar_t_out>() : NULL, epsilon,
                              grad_input->DATA_PTR<scalar_t_in>(),
                              gamma != NULL ? grad_gamma->DATA_PTR<scalar_t_out>() : NULL,
                              beta != NULL ? grad_beta->DATA_PTR<scalar_t_out>() : NULL);)
}