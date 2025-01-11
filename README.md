# Protenix: Protein + X

A trainable PyTorch reproduction of [AlphaFold 3](https://www.nature.com/articles/s41586-024-07487-w).

For more information on the model's performance and capabilities, see our technical report ([biorxiv](https://www.biorxiv.org/content/10.1101/2025.01.08.631967v) | [pdf](Protenix_Technical_Report.pdf)). 

You can follow our [twitter](https://x.com/ai4s_protenix) or join the conversation in the [discord server](https://discord.gg/8ZMWy89aMf).

![Protenix predictions](assets/protenix_predictions.gif)

## ⚡ Try it online
- [Web server link](http://101.126.11.40:8000/) 

## Installation

### Run with PyPI (recommended):

```bash
pip3 install protenix
```
### Run with Docker:

If you're interested in model training, we recommand to [<u> run with docker</u>](docs/docker_installation.md).

## Inference

### Command line inference

If you set up `Protenix` by `pip`, you can run the following command to do model inference:

```bash
# run with example.json, which contains precomputed msa dir.
protenix predict --input examples/example.json --out_dir  ./output --seeds 101

# run with multiple json files, the default seed is 101.
protenix predict --input ./jsons_dir/ --out_dir  ./output

# if the json do not contain precomputed msa dir, 
# add --use_msa_server to search msa and then predict.
# if mutiple seeds are provided, split them by comma.
protenix predict --input examples/example_without_msa.json --out_dir ./output --seeds 101,102 --use_msa_server
```

**Detailed information on the format of the input JSON file and the output files can be found in [<u> input and output documentation </u>](docs/infer_json_format.md)**.

Alternatively you can run inference by:

Note: by default, we do not use layernorm and EvoformerAttention kernels for simple configuration, if you want to speed up inference, see [<u> setting up kernels documentation </u>](docs/kernels.md).

```bash
bash inference_demo.sh
```

Arguments in this scripts are explained as follows:

* `input_json_path`: path to a JSON file that fully describes the input.
* `dump_dir`: path to a directory where the results of the inference will be saved. 
* `dtype`: data type used in inference. Valid options include `"bf16"` and `"fp32"`. 
* `use_msa`: whether to use the MSA feature, the default is true.


### Convert PDB/CIF file to json

If your input is pdb or cif file, you can convert it to json file for inference.
```bash
# run with pdb/cif file, and convert it to json file for inference.
protenix tojson --input examples/7pzb.pdb --out_dir ./output
```


### MSA search
We also provide an independent MSA search function, you can do msa search from json file or fasta file.
```bash
# run msa search with json file, it will write precomputed msa dir info to a new json file.
protenix msa --input examples/example_without_msa.json --out_dir ./output

# run msa search with fasta file which only contains protein.
protenix msa --input examples/prot.fasta --out_dir ./output
```

### Run with PyMol

If you want to run Protenix inference with `PyMol`, please refer to [PyMOLfold](https://github.com/colbyford/PyMOLfold).

## Training
If you're interested in model training, see [<u> training documentation </u>](docs/training.md).

## Performance
See the [<u>performance documentation</u>](docs/model_performance.md) for memory and time consumption in training and inference.

## Citing This Work

If you use this code or the model in your research, please cite the following paper:

```
@article{chen2025protenix,
  title={Protenix - Advancing Structure Prediction Through a Comprehensive AlphaFold3 Reproduction},
  author={Chen, Xinshi and Zhang, Yuxuan and Lu, Chan and Ma, Wenzhi and Guan, Jiaqi and Gong, Chengyue and Yang, Jincai and Zhang, Hanyu and Zhang, Ke and Wu, Shenghao and Zhou, Kuangqi and Yang, Yanping and Liu, Zhenyu and Wang, Lan and Shi, Bo and Shi, Shaochen and Xiao, Wenzhi},
  year={2025},
  doi = {10.1101/2025.01.08.631967},
  journal = {bioRxiv}
}
```


## Acknowledgements

Implementation of the layernorm operators referred to [OneFlow](https://github.com/Oneflow-Inc/oneflow) and [FastFold](https://github.com/hpcaitech/FastFold). We used [OpenFold](https://github.com/aqlaboratory/openfold) for some [module](protenix/openfold_local/) implementations, except the [`LayerNorm`](protenix/model/layer_norm/).


## Contribution

Please check [Contributing](CONTRIBUTING.md) for more details. If you encounter problems using Protenix, feel free to create an issue! We also welcome pull requests from the community.

## Code of Conduct

Please check [Code of Conduct](CODE_OF_CONDUCT.md) for more details.

## Security

If you discover a potential security issue in this project, or think you may
have discovered a security issue, we ask that you notify Bytedance Security via our [security center](https://security.bytedance.com/src) or [vulnerability reporting email](sec@bytedance.com).

Please do **not** create a public GitHub issue.

## License

The Protenix project, including code and model parameters, is made available under the [Apache 2.0 License](./LICENSE), it is free for both academic research and commercial use.

We welcome inquiries and collaboration opportunities for advanced applications of our model, such as developing new features, fine-tuning for specific use cases, and more. Please feel free to contact us at ai4s-bio@bytedance.com.
