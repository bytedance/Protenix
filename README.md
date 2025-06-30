# Protenix: Protein + X


> 📣📣📣 **We're hiring!** \
> Positions in **_Beijing, China_** and **_Seattle, US_** ! \
> Interested in machine learning, computational chemistry/biology, structural biology, or drug discovery?  \
> 👉 [**Join us »**](#join-us)



<div align="center" style="margin: 20px 0;">
  <span style="margin: 0 10px;">⚡ <a href="https://protenix-server.com">Protenix Web Server</a></span>
  &bull; <span style="margin: 0 10px;">📄 <a href="https://www.biorxiv.org/content/10.1101/2025.01.08.631967v1">Technical Report</a></span>
</div>

<div align="center">

[![Twitter](https://img.shields.io/badge/Twitter-Follow-blue?logo=x)](https://x.com/ai4s_protenix)
[![Slack](https://img.shields.io/badge/Slack-Join-yellow?logo=slack)](https://join.slack.com/t/protenixworkspace/shared_invite/zt-36j4kx1cy-GyQMWLDrMO4Wd0fjGxtxug)
[![Wechat](https://img.shields.io/badge/Wechat-Join-brightgreen?logo=wechat)](https://github.com/bytedance/Protenix/issues/52)
[![Email](https://img.shields.io/badge/Email-Contact-lightgrey?logo=gmail)](#contact-us)
</div>

We’re excited to introduce **Protenix** — a trainable, open-source PyTorch reproduction of [AlphaFold 3](https://www.nature.com/articles/s41586-024-07487-w).

Protenix is built for high-accuracy structure prediction. It serves as an initial step in our journey toward advancing accessible and extensible research tools for the computational biology community.



![Protenix predictions](assets/protenix_predictions.gif)

## 🌟 Related Projects
- **[Protenix-Dock](https://github.com/bytedance/Protenix-Dock)**: Our implementation of a classical protein-ligand docking framework that leverages empirical scoring functions. Without using deep neural networks, Protenix-Dock delivers competitive performance in rigid docking tasks.

## Updates
### 🎉 Model Update

- 2025-05-30: **Protenix-v0.5.0** is now available! You may try Protenix-v0.5.0 by accessing the [server](https://protenix-server.com), or upgrade to the latest version using pip.

### 🔥 Feature Update
- 2025-01-16: The preview version of [constraint feature](./README.md#early-access-to-new-constraint-feature) is released to branch [`constraint_esm`](https://github.com/bytedance/Protenix/tree/constraint_esm).
- 2025-01-16: The [training data pipeline](./docs/prepare_training_data.md) is released.
- 2025-01-16: The [MSA pipeline](./docs/msa_pipeline.md) is released.
- 2025-01-16: Use [local colabfold_search](./docs/colabfold_compatible_msa.md) to generate protenix-compatible MSA.

### 📊 Benchmark
We benchmarked the performance of Protenix-v0.5.0 against [Boltz-1](https://github.com/jwohlwend/boltz/releases/tag/v0.4.1) and [Chai-1](https://github.com/chaidiscovery/chai-lab/releases/tag/v0.6.1) across multiple datasets, including [PoseBusters v2](https://arxiv.org/abs/2308.05777), [AF3 Nucleic Acid Complexes](https://www.nature.com/articles/s41586-024-07487-w), [AF3 Antibody Set](https://github.com/google-deepmind/alphafold3/blob/20ad0a21eb49febcaad4a6f5d71aa6b701512e5b/docs/metadata_antibody_antigen.csv), and our curated Recent PDB set.
<!-- 1️⃣ [PoseBusters v2](https://arxiv.org/abs/2308.05777)\
2️⃣ [AF3 Nucleic Acid Complexes](https://www.nature.com/articles/s41586-024-07487-w)\
3️⃣ [AF3 Antibody Set](https://github.com/google-deepmind/alphafold3/blob/20ad0a21eb49febcaad4a6f5d71aa6b701512e5b/docs/metadata_antibody_antigen.csv)\
4️⃣ Our curated Recent PDB set -->

Protenix-v0.5.0 was trained using a PDB cut-off date of September 30, 2021. For the comparative analysis, we adhered to AF3’s inference protocol, generating 25 predictions by employing 5 model seeds, with each seed yielding 5 diffusion samples. The predictions were subsequently ranked based on their respective ranking scores.


![V0.5.0 model Metrics](assets/v0.5.0_metrics.png)

We will soon release the benchmarking toolkit, including the evaluation datasets, data curation pipeline, and metric calculators, to support transparent and reproducible benchmarking.


## 🛠 Installation

### PyPI

```bash
pip3 install protenix
```

For development on a CPU-only machine, it is convenient to install with the `--cpu` flag in editable mode:
```
python3 setup.py develop --cpu
```

### Docker (Recommended for Training)

Check the detailed guide: [<u> Docker Installation</u>](docs/docker_installation.md).


## 🚀 Inference

### Expected Input & Output Format
For details on the input JSON format and expected outputs, please refer to the [Input/Output Documentation](docs/infer_json_format.md).


### Prepare Inputs

#### Convert PDB/CIF File to Input JSON

If your input is a `.pdb` or `.cif` file, you can convert it into a JSON file for inference.


```bash
# ensure `release_data/ccd_cache/components.cif` or run:
python scripts/gen_ccd_cache.py -c release_data/ccd_cache/ -n [num_cpu]

# for PDB
# download pdb file
wget https://files.rcsb.org/download/7pzb.pdb
# run with pdb/cif file, and convert it to json file for inference.
protenix tojson --input examples/7pzb.pdb --out_dir ./output

# for CIF (same process)
# download cif file
wget https://files.rcsb.org/download/7pzb.cif
# run with pdb/cif file, and convert it to json file for inference.
protenix tojson --input examples/7pzb.cif --out_dir ./output
```


#### (Optional) Prepare MSA Files

We provide an independent MSA search utility. You can run it using either a JSON file or a protein FASTA file.
```bash
# run msa search with json file, it will write precomputed msa dir info to a new json file.
protenix msa --input examples/example_without_msa.json --out_dir ./output

# run msa search with fasta file which only contains protein.
protenix msa --input examples/prot.fasta --out_dir ./output
```

### Inference via Command Line

If you installed `Protenix` via `pip`, you can run the following command to perform model inference:


```bash
# the default n_cycle/n_step/n_samples is 10/200/5 respectively, you can modify it by passing --cycle x1 --step x2 --sample x3

# run with example.json, which contains precomputed msa dir.
protenix predict --input examples/example.json --out_dir  ./output --seeds 101

# run with multiple json files, the default seed is 101.
protenix predict --input ./jsons_dir/ --out_dir  ./output

# if the json do not contain precomputed msa dir,
# add --use_msa_server to search msa and then predict.
# if mutiple seeds are provided, split them by comma.
protenix predict --input examples/example_without_msa.json --out_dir ./output --seeds 101,102 --use_msa_server
```

### Inference via Bash Script
Alternatively you can run inference by:
Alternatively, run inference via script:

```bash
bash inference_demo.sh
```

The script accepts the following arguments:
* `input_json_path`: Path to a JSON file that fully specifies the input structure.
* `dump_dir`: Directory where inference results will be saved.
* `dtype`: Data type used during inference. Supported options: `bf16` and `fp32`.
* `use_msa`: Whether to enable MSA features (default: true).
* `use_esm`: Whether to enable ESM features (default: false).


> **Note**: By default, layernorm and EvoformerAttention kernels are disabled for simplicity.
> To enable them and speed up inference, see the [**Kernels Setup Guide**](docs/kernels.md).


## 🧬 Training

Refer to the [Training Documentation](docs/training.md) for setup and details.

## 📌 Constraint Feature

Protenix now allows users to specify ***contacts***, enabling the model to leverage additional inter-chain information as constraint guidance! We benchmarked this feature on the PoseBusters dataset and a curated protein-antibody interface subset.  Results show that Protenix can generate significantly more accurate structures when guided by constraints. You can try it out via the [`constraint_esm`](https://github.com/bytedance/Protenix/tree/constraint_esm) branch.

![Constraint Metrics](assets/constraint_metrics.png)

> **Tips:** Our online service already supports constraint inputs — no local setup required!
However, for local command-line usage, be sure to check out the [`constraint_esm`](https://github.com/bytedance/Protenix/tree/constraint_esm) branch, as this feature is not yet included in the main branch.


## Training and Inference Cost

For details on memory usage and runtime during training and inference, refer to the [Training & Inference Cost Documentation](docs/model_train_inference_cost.md).


## Citing Protenix

If you use Protenix in your research, please cite the following:

```
@article{bytedance2025protenix,
  title={Protenix - Advancing Structure Prediction Through a Comprehensive AlphaFold3 Reproduction},
  author={ByteDance AML AI4Science Team and Chen, Xinshi and Zhang, Yuxuan and Lu, Chan and Ma, Wenzhi and Guan, Jiaqi and Gong, Chengyue and Yang, Jincai and Zhang, Hanyu and Zhang, Ke and Wu, Shenghao and Zhou, Kuangqi and Yang, Yanping and Liu, Zhenyu and Wang, Lan and Shi, Bo and Shi, Shaochen and Xiao, Wenzhi},
  year={2025},
  journal={bioRxiv},
  publisher={Cold Spring Harbor Laboratory},
  doi={10.1101/2025.01.08.631967},
  URL={https://www.biorxiv.org/content/early/2025/01/11/2025.01.08.631967},
  elocation-id={2025.01.08.631967},
  eprint={https://www.biorxiv.org/content/early/2025/01/11/2025.01.08.631967.full.pdf},
}
```

### 📚 Citing Related Work
Protenix is built upon and inspired by several influential projects. If you use Protenix in your research, we also encourage citing the following foundational works where appropriate:
```
@article{abramson2024accurate,
  title={Accurate structure prediction of biomolecular interactions with AlphaFold 3},
  author={Abramson, Josh and Adler, Jonas and Dunger, Jack and Evans, Richard and Green, Tim and Pritzel, Alexander and Ronneberger, Olaf and Willmore, Lindsay and Ballard, Andrew J and Bambrick, Joshua and others},
  journal={Nature},
  volume={630},
  number={8016},
  pages={493--500},
  year={2024},
  publisher={Nature Publishing Group UK London}
}
@article{ahdritz2024openfold,
  title={OpenFold: Retraining AlphaFold2 yields new insights into its learning mechanisms and capacity for generalization},
  author={Ahdritz, Gustaf and Bouatta, Nazim and Floristean, Christina and Kadyan, Sachin and Xia, Qinghui and Gerecke, William and O’Donnell, Timothy J and Berenberg, Daniel and Fisk, Ian and Zanichelli, Niccol{\`o} and others},
  journal={Nature Methods},
  volume={21},
  number={8},
  pages={1514--1524},
  year={2024},
  publisher={Nature Publishing Group US New York}
}
@article{mirdita2022colabfold,
  title={ColabFold: making protein folding accessible to all},
  author={Mirdita, Milot and Sch{\"u}tze, Konstantin and Moriwaki, Yoshitaka and Heo, Lim and Ovchinnikov, Sergey and Steinegger, Martin},
  journal={Nature methods},
  volume={19},
  number={6},
  pages={679--682},
  year={2022},
  publisher={Nature Publishing Group US New York}
}
```

## Contributing to Protenix

We welcome contributions from the community to help improve Protenix!

📄 Check out the [Contributing Guide](CONTRIBUTING.md) to get started.

✅ Code Quality: 
We use `pre-commit` hooks to ensure consistency and code quality. Please install them before making commits:

```bash
pip install pre-commit
pre-commit install
```

🐞 Found a bug or have a feature request? [Open an issue](https://github.com/bytedance/Protenix/issues).



## Acknowledgements


The implementation of LayerNorm operators refers to both [OneFlow](https://github.com/Oneflow-Inc/oneflow) and [FastFold](https://github.com/hpcaitech/FastFold).
We also adopted several [module](protenix/openfold_local/) implementations from [OpenFold](https://github.com/aqlaboratory/openfold), except for [`LayerNorm`](protenix/model/layer_norm/), which is implemented independently.


## Code of Conduct

We are committed to fostering a welcoming and inclusive environment.
Please review our [Code of Conduct](CODE_OF_CONDUCT.md) for guidelines on how to participate respectfully.


## Security

If you discover a potential security issue in this project, or think you may
have discovered a security issue, we ask that you notify Bytedance Security via our [security center](https://security.bytedance.com/src) or [vulnerability reporting email](sec@bytedance.com).

Please do **not** create a public GitHub issue.

## License

The Protenix project including both code and model parameters is released under the [Apache 2.0 License](./LICENSE). It is free for both academic research and commercial use.

## Contact Us

We welcome inquiries and collaboration opportunities for advanced applications of our model, such as developing new features, fine-tuning for specific use cases, and more. Please feel free to contact us at ai4s-bio@bytedance.com.

## Join Us

We're expanding the **Protenix team** at ByteDance Seed-AI for Science! \
We’re looking for talented individuals in **machine learning** and **computational biology/chemistry**. Opportunities are available in both **Beijing** and **Seattle**, across internships, new grad roles, and experienced full-time positions. \
*“Computational Biology/Chemistry” covers structural biology, computational biology, computational chemistry, drug discovery, and more.*


### 📍 Beijing, China
| Type       | Expertise                          | Apply Link |
|------------|------------------------------------|------------|
| Full-Time  | Computational Biology / Chemistry       | [Experienced & New Grad](https://jobs.bytedance.com/society/position/detail/7505998274429421842) |
| Full-Time  | Machine Learning                   | [Experienced & New Grad](https://jobs.bytedance.com/society/position/detail/7505999453133015314) |
| Internship | Computational Biology / Chemistry       | [Internship](https://jobs.bytedance.com/campus/position/7509005072577546504/detail) |
| Internship | Machine Learning                   | [Internship](https://jobs.bytedance.com/campus/position/7509005074018961672/detail) |


### 📍 Seattle, US

| Type       | Expertise                          | Apply Link |
|------------|------------------------------------|------------|
| Full-Time  | Computational Biology / Chemistry       | [Experienced](https://jobs.bytedance.com/en/position/7270666468370614585/detail), [New Grad](https://jobs.bytedance.com/en/position/7515465250054211847/detail) |
| Full-Time  | Machine Learning                   | [Experienced](https://jobs.bytedance.com/en/position/7270665658072926521/detail), [New Grad](https://jobs.bytedance.com/en/position/7515908698011601159/detail) |
| Internship | Computational Biology / Chemistry       | Internship (opening ~August) |
| Internship | Machine Learning                   | Internship (opening ~August) |
