# Protenix: Protein + X

A trainable PyTorch reproduction of [AlphaFold 3](https://www.nature.com/articles/s41586-024-07487-w).

For more information on the model's performance and capabilities, see our [technical report](Protenix_Technical_Report.pdf). 

You can follow our [twitter](https://x.com/ai4s_protenix) or join the conversation in the [discord server](https://discord.gg/Utgk4Ykw).

![Protenix predictions](assets/protenix_predictions.gif)

## ⚡ Try it online
- [Web server link](http://101.126.11.40:8000/) 

## Installation

### Run with PyPI (recommended):

```bash
    # maybe you need to update libxrender1 and libxext6 firstly, run as following for Debian:
    # apt-get update && apt-get install -y libxrender1 libxext6
    pip3 install protenix
```
### Run with Docker:

If you're interested in model training, we recommand to [<u> run with docker</u>](docs/docker_installation.md).

## Inference

### Command line inference

If you set up `Protenix` by `pip`, you can run the following command to do model inference:

```bash
# run with one json file
protenix predict --input examples/example.json --out_dir  ./output

# or run with multiple json files
protenix predict --input ./jsons_dir/ --out_dir  ./output

```

**Detailed information on the format of the input JSON file and the output files can be found in [<u> input and output documentation </u>](docs/infer_json_format.md)**.

Alternatively you can run inference by:

```bash
bash inference_demo.sh
```

Arguments in this scripts are explained as follows:
* `load_checkpoint_path`: path to the model checkpoints.
* `input_json_path`: path to a JSON file that fully describes the input.
* `dump_dir`: path to a directory where the results of the inference will be saved. 
* `dtype`: data type used in inference. Valid options include `"bf16"` and `"fp32"`. 
* `use_msa`: whether to use the MSA feature, the default is true. If you want to disable the MSA feature, add `--use_msa false` to the [inference_demo.sh](inference_demo.sh) script.


Note: by default, we do not use layernorm and EvoformerAttention kernels for simple configuration, if you want to speed up inference, see [<u> setting up kernels documentation </u>](docs/kernels.md).

### Notebook demo
You can use [notebooks/protenix_inference.ipynb](notebooks/protenix_inference.ipynb)  to run the model inference.

## Training
If you're interested in model training, see [<u> training documentation </u>](docs/training.md).

## Performance
See the [<u>performance documentation</u>](docs/model_performance.md) for memory and time consumption in training and inference.

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

This project, including code and model parameters are made available under the terms of the Creative Commons Attribution-NonCommercial 4.0 International License. You can find details at: https://creativecommons.org/licenses/by-nc/4.0/

For commercial use, please reach out to us at ai4s-bio@bytedance.com for the commercial license. We welcome all types of collaborations.
