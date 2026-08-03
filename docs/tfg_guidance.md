# Training-Free Guidance for User Contact Constraints

Training-Free Guidance (TFG) applies differentiable energy terms during
diffusion sampling. It is a sampling-time guidance mechanism, not a learned
conditioning module.

## Constraint Support by Model

- `protenix_base_constraint_v0.5.0` uses the trained constraint embedder path
  for JSON `constraint` inputs.
- `protenix_base_default_v1.0.0`, `protenix_base_20250630_v1.0.0`, and
  `protenix-v2` do not have trained constraint embedders enabled.
- For v1/v2 models, JSON `constraint.contact` entries can be used as TFG
  distance restraints by enabling TFG.
- JSON `constraint.pocket` entries are not converted to TFG restraints yet for
  v1/v2 models. Use the v0.5 constraint model for pocket constraints.

TFG restraints are soft. They can bias sampling toward satisfying the specified
distance bounds, but they do not guarantee satisfaction and can trade off with
model confidence or geometry if weighted too strongly.

## Contact Example

Use the existing JSON `constraint.contact` format:

```json
[
  {
    "name": "example_tfg_contact",
    "sequences": [
      {
        "proteinChain": {
          "sequence": "AAAA",
          "count": 1
        }
      },
      {
        "proteinChain": {
          "sequence": "GGGG",
          "count": 1
        }
      }
    ],
    "constraint": {
      "contact": [
        {
          "entity1": 1,
          "copy1": 1,
          "position1": 2,
          "entity2": 2,
          "copy2": 1,
          "position2": 3,
          "max_distance": 8.0,
          "min_distance": 0.0
        }
      ]
    }
  }
]
```

Run with the Click CLI:

```bash
protenix pred \
  -i examples/example_tfg_contact_v1.json \
  -o ./output_tfg_contact \
  -n protenix_base_default_v1.0.0 \
  --use_tfg_guidance true \
  --tfg_constraint_weight 1.0
```

Advanced users can tune through nested config overrides in `runner/inference.py`:

```bash
python3 runner/inference.py \
  --model_name protenix_base_default_v1.0.0 \
  --input_json_path examples/example_tfg_contact_v1.json \
  --dump_dir ./output_tfg_contact \
  --seeds 101 \
  --sample_diffusion.N_sample 1 \
  --sample_diffusion.N_step 200 \
  --sample_diffusion.guidance.enable true \
  --sample_diffusion.guidance.terms.UserDistanceRestraintPotential.weight 1.0
```

## Atom-Level Contacts

For exact ligand or residue atom control, specify both atom names:

```json
{
  "entity1": 1,
  "copy1": 1,
  "position1": 10,
  "atom1": "CA",
  "entity2": 2,
  "copy2": 1,
  "position2": 1,
  "atom2": "C5",
  "max_distance": 6.0,
  "min_distance": 3.0
}
```

When atom names are omitted, Protenix resolves the contact to deterministic
token centre atoms. For ligands, atom-level contacts are recommended because
ligands are represented by atom-level tokens.

## Tuning

- Start with `--sample_diffusion.N_sample 1` while tuning.
- Increase `--tfg_constraint_weight` if restraints are consistently violated.
- Decrease the weight if structures become distorted or confidence drops.
- TFG adds work inside diffusion sampling, so guided runs are slower than
  unguided runs.

## Troubleshooting

- If v1/v2 constraints are present but `--use_tfg_guidance` is false, they are
  parsed but do not affect sampling.
- If `constraint.pocket` is present with a v1/v2 model, it is logged as not yet
  converted to TFG restraints.
- If a contact cannot be resolved to input atoms, it is skipped and logged.
- If duplicate contact entries give incompatible bounds for the same atom pair,
  inference raises a validation error.
