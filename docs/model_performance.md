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

The model will be infered in BF16 Mixed precision, by **default**, the `SampleDiffusion`,`ConfidenceHead` part will still be infered in FP32 precision. if you want to infer the model in **full BF16** Mixed precision, pass the following arguments to the [inference_demo.sh](../inference_demo.sh):

  ```bash
  --skip_amp.confidence_head false \
  --skip_amp.sample_diffusion false \
  ```

Below are reference examples of cuda memory usage (G).

| Ntoken | Natom | Default | Full BF16 Mixed |
|--------|-------|-------|------------------|
| 500    | 10000 | 5.6   | 5.1  |
| 900    | 18000 | 10.7  | 8.1  |
| 1100   | 22000 | 14.4  | 12.1 |
| 1300   | 26000 | 22.1  | 15.0 |
| 1500   | 30000 | 24.8  | 19.2 |
| 1700   | 34000 | 30.3  | 24.1 |
| 1900   | 38000 | 37.1  | 29.3 |
| 2000   | 20000 | 48.0  | 26.5 |
| 2500   | 25000 | 52.2  | 34.8 |
| 3000   | 30000 | 73.8  | 44.0 |
| 3500   | 35000 | 67.6  | 38.2 |
| 4000   | 40000 | 66.4  | 48.4 |
| 4500   | 45000 | 77.0  | 59.2 |
| 5000   | 50000 | OOM   | 72.8 |