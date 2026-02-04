# Protenix v0.5.0 Series Benchmark & Model Features

This document provides benchmark results and descriptions of key features for the Protenix v0.5.0 series models.

## 0.5.0 Series Models
The v0.5.0 series includes the following model variants:
- `protenix_base_default_v0.5.0`: Standard base model.
- `protenix_base_constraint_v0.5.0`: Base model with constraint support.
- `protenix_mini_esm_v0.5.0`: Lightweight model using ESM features.
- `protenix_mini_ism_v0.5.0`: Lightweight model using ISM features.
- `protenix_mini_default_v0.5.0`: Lightweight model without extra ESM/ISM features.
- `protenix_tiny_default_v0.5.0`: Smallest model variant for maximum efficiency.

---

## 📊 Benchmark
We benchmarked the performance of Protenix-v0.5.0 against [Boltz-1](https://github.com/jwohlwend/boltz/releases/tag/v0.4.1) and [Chai-1](https://github.com/chaidiscovery/chai-lab/releases/tag/v0.6.1) across multiple datasets, including [PoseBusters v2](https://arxiv.org/abs/2308.05777), [AF3 Nucleic Acid Complexes](https://www.nature.com/articles/s41586-024-07487-w), [AF3 Antibody Set](https://github.com/google-deepmind/alphafold3/blob/20ad0a21eb49febcaad4a6f5d71aa6b701512e5b/docs/metadata_antibody_antigen.csv), and our curated Recent PDB set.

Protenix-v0.5.0 was trained using a PDB cut-off date of September 30, 2021. For the comparative analysis, we adhered to AF3’s inference protocol, generating 25 predictions by employing 5 model seeds, with each seed yielding 5 diffusion samples. The predictions were subsequently ranked based on their respective ranking scores.

![V0.5.0 model Metrics](../assets/protenix_base_default_v0.5.0_metrics.png)

We will soon release the benchmarking toolkit, including the evaluation datasets, data curation pipeline, and metric calculators, to support transparent and reproducible benchmarking.

---

## Model Features

### 📌 Constraint
Protenix supports specifying ***contacts*** (at both residue and atom levels) and ***pocket constraints*** as extra guidance. Our benchmark results demonstrate that constraint-guided predictions are significantly more accurate. See our [doc](./infer_json_format.md#constraint) for input format details.

![Constraint Metrics](../assets/protenix_base_constraint_v0.5.0_metrics.png)

### 📌 Mini-Models
We introduce Protenix-Mini, a lightweight variant of Protenix that uses reduced network blocks and few ODE steps (even as few as one or two steps) to enable efficient prediction of biomolecular complex structures. Experimental results show that Protenix-Mini achieves a favorable balance between efficiency and accuracy, with only a marginal 1–5% drop in evaluation metrics such as interface LDDT, complex LDDT, and ligand RMSD success rate. Protenix-Mini enables accurate structure prediction in high-throughput and resource-limited scenarios, making it well-suited for practical applications at scale. The following comparisons were performed on a subset of the RecentPDB dataset comprising sequences with fewer than 768 tokens.

![Mini/Tiny Metrics](../assets/mini_tiny_0.5.0_performance.png)
