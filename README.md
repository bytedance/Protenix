# Protenix: Protein + X


> 📣📣📣 **We're hiring!** \
> Positions in **_Beijing_** 🇨🇳 and **_Seattle_** 🇺🇸 \
> Interested in machine learning, computational chemistry/biology, structural biology, or drug discovery?  Join us to build cutting-edge AI for biology!
\
> 👉 [**Join us »**](#join-us)



<div align="center" style="margin: 20px 0;">
  <span style="margin: 0 10px;">⚡ <a href="https://protenix-server.com">Protenix Web Server</a></span>
  &bull; <span style="margin: 0 10px;">📄 <a href="https://www.biorxiv.org/content/10.1101/2025.01.08.631967v1">Technical Report</a></span>
</div>

<div align="center">

[![Twitter](https://img.shields.io/badge/Twitter-Follow-blue?logo=x)](https://x.com/ai4s_protenix)
[![Slack](https://img.shields.io/badge/Slack-Join-yellow?logo=slack)](https://join.slack.com/t/protenixworkspace/shared_invite/zt-3drypwagk-zRnDF2VtOQhpWJqMrIveMw)
[![Wechat](https://img.shields.io/badge/Wechat-Join-brightgreen?logo=wechat)](https://github.com/bytedance/Protenix/issues/52)
[![Email](https://img.shields.io/badge/Email-Contact-lightgrey?logo=gmail)](#contact-us)
</div>

We’re excited to introduce **Protenix** — Toward High-Accuracy Open-Source Biomolecular Structure Prediction.

Protenix is built for high-accuracy structure prediction. It serves as an initial step in our journey toward advancing accessible and extensible research tools for the computational biology community.

<img src="assets/protenix_predictions.gif" style="width: 100%; height: auto;" alt="Protenix predictions">

## 🌟 Related Projects
- **[PXDesign](https://protenix.github.io/pxdesign/)** is a model suite for de novo protein-binder design built on the Protenix foundation model. PXDesign achieves 20–73% experimental success rates across multiple targets — 2–6× higher than prior SOTA methods such as AlphaProteo and RFdiffusion. The framework is freely accessible via the Protenix Server.

- **[PXMeter](https://github.com/bytedance/PXMeter/)** is an open-source toolkit designed for reproducible evaluation of structure prediction models, released with high-quality benchmark dataset that has been manually reviewed to remove experimental artifacts and non-biological interactions. The associated study presents an in-depth comparative analysis of state-of-the-art models, drawing insights from extensive metric data and detailed case studies. The evaluation of Protenix is based on PXMeter.

- **[Protenix-Dock](https://github.com/bytedance/Protenix-Dock)**: Our implementation of a classical protein-ligand docking framework that leverages empirical scoring functions. Without using deep neural networks, Protenix-Dock delivers competitive performance in rigid docking tasks.

## 🎉 Latest Updates
- **2026-02-05: Protenix-v1 Released** 💪 [[Technical Report](docs/PTX_V1_Technical_Report_202602042356.pdf)]
  - Supported Template/RNA MSA features and improved training dynamics, along with further Inference-time model performance enhancements.
- **2025-11-05: Protenix-v0.7.0 Released** 🚀
  - Introduced advanced diffusion inference optimizations: Shared variable caching, efficient kernel fusion, and TF32 acceleration. See our [performance analysis](./assets/inference_time_vs_ntoken.png).
- **2025-07-17: Protenix-Mini & Constraint Features**
  - Released lightweight model variants ([Protenix-Mini](https://arxiv.org/abs/2507.11839)) that drastically reduce inference costs with minimal accuracy loss.
  - Added support for [atom-level contact and pocket constraints](docs/infer_json_format.md#constraint), enhancing prediction accuracy through physical priors.
- **2025-01-16: Pipeline Enhancements**
  - Open-sourced the full [training data pipeline](./docs/prepare_training_data.md) and [MSA pipeline](./docs/msa_template_pipeline.md).
  - Integrated local [ColabFold-compatible search](./docs/colabfold_compatible_msa.md) for streamlined MSA generation.


## 🚀 Getting Started

### 🛠 Quick Installation

```bash
pip install protenix
```

### 🧬 Quick Prediction

```bash
# Predict structure using a JSON input
protenix pred -i examples/input.json -o ./output -n protenix_base_default_v1.0.0
```

For detailed instructions on installation, data preprocessing, inference, and training, please refer to the [Training and Inference Instructions](docs/training_inference_instructions.md) and [Supported Models](docs/supported_models.md). We recommend users refer to [inference_demo.sh](inference_demo.sh) for detailed inference methods and input explanations.


### 📊 Benchmark

Protenix-v1, the first fully open-source model that outperforms AlphaFold3 across diverse benchmark sets while adhering to the same training data cutoff, model scale, and inference budget as AlphaFold3. For challenging targets, such as antigen-antibody complexes, the prediction accuracy of Protenix-v1 can be further enhanced through inference-time scaling – increasing the sampling budget from several to hundreds of candidates leads to consistent log-linear gains.

<img src="./assets/protenix_base_default_v1.0.0_metrics.png" style="width: 100%; height: auto;" alt="protenix-v1 model Metrics">

<img src="./assets/protenix_base_default_v1.0.0_metrics2.png" style="width: 100%; height: auto;" alt="protenix-v1 model Metrics 2">

For detailed benchmark metrics on each dataset, please refer to [docs/model_1.0.0_benchmark.md](docs/model_1.0.0_benchmark.md).

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

We're expanding the **Protenix team** at ByteDance Seed-AI for Science! We’re looking for talented individuals in machine learning and computational biology/chemistry (*“Computational Biology/Chemistry” covers structural biology, computational biology, computational chemistry, drug discovery, and more*). Opportunities are available in both **Beijing** and **Seattle**, across internships, new grad roles, and experienced full-time positions. 

Outstanding applicants will be considered for **ByteDance’s Top Seed Talent Program** — with enhanced support.


### 📍 Beijing, China
| Type       | Expertise                          | Apply Link |
|------------|------------------------------------|------------|
| Full-Time  | Protein Design Scientist       | [Experienced](https://jobs.bytedance.com/society/position/detail/7550992796392982792) |
| Full-Time  | Computational Biology / Chemistry       | [Experienced](https://jobs.bytedance.com/society/position/detail/7505998274429421842), [New Grad](https://job.toutiao.com/s/HGwWBs1UGR4) |
| Full-Time  | Machine Learning                   | [Experienced](https://jobs.bytedance.com/society/position/detail/7505999453133015314), [New Grad](https://job.toutiao.com/s/upy82CljXlY) |
| Internship | Computational Biology / Chemistry       | [Internship](https://job.toutiao.com/s/Wr3yig1Wet4) |
| Internship | Machine Learning                   | [Internship](https://job.toutiao.com/s/w2GQQDfQUkc) |


### 📍 Seattle, US

| Type       | Expertise                          | Apply Link |
|------------|------------------------------------|------------|
| Full-Time  | Computational Biology / Chemistry       | [Experienced](https://jobs.bytedance.com/en/position/7270666468370614585/detail), [New Grad](https://job.toutiao.com/s/iH00nSEvrFo) |
| Full-Time  | Machine Learning                   | [Experienced](https://jobs.bytedance.com/en/position/7270665658072926521/detail), [New Grad](https://job.toutiao.com/s/dmU_fbEHGOw) |
| Internship | Computational Biology / Chemistry       | [Internship](https://job.toutiao.com/s/aiCZz0kJexs) |
| Internship | Machine Learning                   | [Internship](https://job.toutiao.com/s/DiGnn5l1QpQ) |
