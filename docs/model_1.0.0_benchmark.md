# Protenix v1.0.0 Series Benchmarks

This document provides detailed benchmark performance results and access to the evaluation datasets used for Protenix v1.0.0.


## 📂 Benchmark Data & Results

The complete benchmark results are available for download. Each package contains detailed inference metrics for individual predictions as well as aggregated performance reports.

| Dataset | Size | Description | Download |
| :--- | :--- | :--- | :--- |
| `PXM-2024` | ~32 MB | Benchmark results for structures deposited in 2024. | [Link](https://protenix.tos-cn-beijing.volces.com/benchmark/v1/PXM-2024.tar.gz) |
| `PXM-2025` | ~28 MB | Benchmark results for structures deposited in 2025. | [Link](https://protenix.tos-cn-beijing.volces.com/benchmark/v1/PXM-2025.tar.gz) |
| `PXM-2025-H2` | ~18 MB | Benchmark results for structures deposited in 2025 (second half). | [Link](https://protenix.tos-cn-beijing.volces.com/benchmark/v1/PXM-2025-H2.tar.gz) |
| `PXM-22to25-Ab-Ag` | ~18 MB | Focused benchmark on Antibody-Antigen complexes (2022-2025). | [Link](https://protenix.tos-cn-beijing.volces.com/benchmark/v1/PXM-22to25-Ab-Ag.tar.gz) |
| `PXM-22to25-Ligand` | ~6 MB | Focused benchmark on Protein-Ligand complexes (2022-2025). | [Link](https://protenix.tos-cn-beijing.volces.com/benchmark/v1/PXM-22to25-Ligand.tar.gz) |
| `PXM-Legacy` | ~15 MB | Historical benchmark sets used in previous development cycles. | [Link](https://protenix.tos-cn-beijing.volces.com/benchmark/v1/PXM-Legacy.tar.gz) |
| `foldbench` | ~11 MB | Evaluation results on the external FoldBench dataset. | [Link](https://protenix.tos-cn-beijing.volces.com/benchmark/v1/foldbench.tar.gz) |
| `runs-n-poses` | ~15 MB | Benchmark focused on protein-ligand and multi-pose evaluation. | [Link](https://protenix.tos-cn-beijing.volces.com/benchmark/v1/runs-n-poses.tar.gz) |

---

## How to Analyze Results

To use these benchmark results:
1. **Download**: Obtain the relevant `.tar.gz` file using the links provided above.
2. **Extract**: Decompress the file using `tar -xzvf <dataset_name>.tar.gz`.
3. **Data Structure**:
   - **Inference Metrics**: Detailed per-prediction metrics (e.g., LDDT, RMSD, confidence scores) provided in CSV format for each sample.
   - **Aggregated Results**: Summary reports providing statistical overviews across the entire dataset for quick performance assessment.
