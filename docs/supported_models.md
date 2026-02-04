# Protenix Supported Models

Protenix provides various pre-trained models to suit different computational resources, inference speeds, and prediction accuracy requirements. This document details the characteristics and configurations of these models.

## Naming Convention

Model names follow the format:
`protenix_{model_size}_{features}_{version}`

- **model_size**: Size of the model, including `base`, `mini` (lightweight), and `tiny` (minimal).
- **features**: Functional characteristics such as `default`, `constraint` (distance constraints), `esm` (includes ESM embeddings), etc. Multiple features are separated by `-`.
- **version**: Version number, e.g., `v0.5.0`, `v1.0.0`.

## Supported Models Summary
| Model Name | ESM | MSA | Constraint | RNA MSA | Template | Params | Training Data Cutoff |
| :--- | :---: | :---: | :---: | :---: | :---: | :---: | :---: |
| `protenix_base_default_v1.0.0` | ❌ | ✅ | ❌ | ✅ | ✅ | 368.48 M | 2021-09-30 |
| `protenix_base_20250630_v1.0.0` * | ❌ | ✅ | ❌ | ✅ | ✅ | 368.48 M | 2025-06-30 |
| `protenix_base_default_v0.5.0` | ❌ | ✅ | ❌ | ❌ | ❌ | 368.09 M | 2021-09-30 |
| `protenix_base_constraint_v0.5.0` | ❌ | ✅ | ✅ | ❌ | ❌ | 368.30 M | 2021-09-30 |
| `protenix_mini_esm_v0.5.0` | ✅ | ✅ | ❌ | ❌ | ❌ | 135.22 M | 2021-09-30 |
| `protenix_mini_ism_v0.5.0` | ✅ | ✅ | ❌ | ❌ | ❌ | 135.22 M | 2021-09-30 |
| `protenix_mini_default_v0.5.0` | ❌ | ✅ | ❌ | ❌ | ❌ | 134.06 M | 2021-09-30 |
| `protenix_tiny_default_v0.5.0` | ❌ | ✅ | ❌ | ❌ | ❌ | 109.50 M | 2021-09-30 |

*\*Note: For practical application scenarios, `protenix_base_20250630_v1.0.0` is trained based on the 2025-06-30 wwPDB cutoff and is also released to the community. For fair benchmarks of model improvements across different versions, please use `protenix_base_default_v1.0.0` (trained on the 2021-09-30 cutoff).*

---

## Model Detailed Descriptions

### 1. Base Models
- **Characteristics**: Full-parameter models with the highest prediction accuracy.
- **Key Configurations**:
    - `N_cycle`: 10 (Number of recycle iterations).
    - `sample_diffusion.N_step`: 200 (Diffusion steps for higher quality sampling).
- **Use Case**: Scientific research requiring maximum precision.

### 2. Mini & Tiny Models
- **Characteristics**: Significant reduction in parameters and faster inference speed.
- **Key Configurations**:
    - `N_cycle`: 4
    - `sample_diffusion.N_step`: 5
- **Difference**: `Mini` has more layers in Pairformer and Transformer modules compared to `Tiny`.
- **Use Case**: High-throughput screening or scenarios with limited computational resources.

### 3. Constraint Model (`protenix_base_constraint_v0.5.0`)
- **Characteristics**: Allows incorporating additional experimental constraints during inference (e.g., Pocket, Contact).
- **Included Features**:
    - `pocket_embedder`: Handles binding pocket information.
    - `contact_embedder`: Handles contact point information.
- **Use Case**: Predictions with available structural priors.

### 4. ESM & ISM Models
- **Characteristics**: Integrates the single-sequence protein language model (ESM2-3B), performing better when MSAs are unavailable.
- **Difference**: `ESM` uses standard ESM2 embeddings, while `ISM` uses specific ISM embeddings.
- **Note**: For efficiency, these models do not use MSA by default.