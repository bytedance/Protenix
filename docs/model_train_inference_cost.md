### Training 
Some settings follow those in the [AlphaFold 3](https://www.nature.com/articles/s41586-024-07487-w) paper, The table below shows the training settings for different fine-tuning stages:

  | Arguments  | Initial training | Fine tuning 1   |  Fine tuning 2  | Fine tuning 3 |
  |-----------------------------------------|--------|---------|-------|-----|
  | `train_crop_size`                       | 384    | 640    | 768    | 768 |
  | `diffusion_batch_size`                  | 48     | 32     | 32     | 32  |
  | `loss.weight.alpha_pae`                 | 0      | 0      | 0      | 1.0 |
  | `loss.weight.alpha_bond`                | 0      | 1.0    | 1.0    | 0   | 
  | `loss.weight.smooth_lddt`               | 1.0    | 0      | 0      | 0   | 
  | `loss.weight.alpha_confidence`          | 1e-4   | 1e-4   | 1e-4   | 1e-4|
  | `loss.weight.alpha_diffusion`           | 4.0    | 4.0    | 4.0    | 0   |
  | `loss.weight.alpha_distogram`           | 0.03   | 0.03   | 0.03   | 0   |
  | `train_confidence_only`                 | False  | False  | False  | True|
  | full BF16-mixed speed(A100, s/step)     | ~12    | ~30    | ~44    | ~13 |
  | full BF16-mixed peak memory (G)         | ~34    | ~35    | ~48    | ~24 |
  
  We recommend carrying out the training on A100-80G or H20/H100 GPUs. If utilizing full BF16-Mixed precision training, the initial training stage can also be performed on A800-40G GPUs. GPUs with smaller memory, such as A30, you'll need to reduce the model size, such as decreasing `model.pairformer.nblocks` and `diffusion_batch_size`.

### Inference

The model will be infered in BF16 Mixed precision, by **default**, the `SampleDiffusion`,`ConfidenceHead` part will still be infered in FP32 precision.

Below are reference examples of cuda memory usage (G).

| Ntoken | Natom | Memory(G) | Inference time(s) |
|--------|-------|-------|------------------|
| 500    | 5000  | 5.2   | 72   |
| 1000   | 10000 | 11.5  | 229  |
| 2000   | 20000 | 42.8  | 933  |
| 3000   | 30000 | 73.9  | 2295 |
| 3500   | 35000 | 69.5  | 3329 |
| 4000   | 40000 | 67.5  | 4483 |

The script in [runner/inference.py](../runner/inference.py) will automatically change the default precision to compute `SampleDiffusion`,`ConfidenceHead` to avoid OOM as follows:
```python
def update_inference_configs(configs: Any, N_token: int):
    # Setting the default inference configs for different N_token and N_atom
    # when N_token is larger than 3000, the default config might OOM even on a
    # A100 80G GPUS,
    if N_token > 3840:
        configs.skip_amp.confidence_head = False
        configs.skip_amp.sample_diffusion = False
    elif N_token > 2560:
        configs.skip_amp.confidence_head = False
        configs.skip_amp.sample_diffusion = True
    else:
        configs.skip_amp.confidence_head = True
        configs.skip_amp.sample_diffusion = True
    return configs
```