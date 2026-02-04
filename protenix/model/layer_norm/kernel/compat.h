// modified from https://github.com/NVIDIA/apex/blob/master/csrc/compat.h
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

#ifndef TORCH_CHECK
#define TORCH_CHECK AT_CHECK
#endif

#ifdef VERSION_GE_1_3
#define DATA_PTR data_ptr
#else
#define DATA_PTR data
#endif
