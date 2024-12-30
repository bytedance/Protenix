### Format of the input JSON file
The JSON file format closely resembles that used by the AlphaFold Server, with a few key differences:

1. There are no restrictions on the types of ligands, ions, and modifications, whereas the AlphaFold Server currently supports only a limited set of specific CCD codes.
2. Users can specify bonds between entities, such as covalent bonds between ligands and polymers.
3. It supports inputting ligands in the form of SMILES strings or molecular structure files.
4. Ligands composed of multiple CCD codes can be treated as a single entity. This feature is useful for representing glycans, for example, "NAG-NAG".
5. The "glycans" field is no longer supported. Glycans can be fully represented by inputting multiple ligands with defined bonding or by providing their SMILES strings.

Here is an overview of the JSON file format:
```json
[ 
  {
    "name": "Test Fold Job Number One",
    "sequences": [...],
    "covalent_bonds": [...]
  }
]
```
The JSON file consists of a list of dictionaries, where each dictionary represents a set of sequences you want to model. 
Even if you are modeling only one set of sequences, the top-level structure should still be a list.

Each dictionary contains the following three keys:
* `name`: A string representing the name of the inference job.
* `sequences`: A list of dictionaries that describe the entities (e.g., proteins, DNA, RNA, small molecules, and ions) involved in the inference.
* `covalent_bonds`: An optional list of dictionaries that define the covalent bonds between atoms from different entities.

Details of `sequences` and `covalent_bonds` are provided below.

#### sequences
There are 5 kinds of supported sequences:
*   `proteinChain` – used for proteins
*   `dnaSequence` – used for DNA (single strand)
*   `rnaSequence` – used for RNA (single strand)
*   `ligand` – used for ligands
*   `ion` – used for ions

##### proteinChain
```json
{
  "proteinChain": {
    "sequence": "PREACHINGS", 
    "count": 1,
    "modifications": [
      {
        "ptmType": "CCD_HY3", 
        "ptmPosition": 1,
      },
      {
        "ptmType": "CCD_P1L",
        "ptmPosition": 5
      }
    ],
    "msa":{
      "precomputed_msa_dir": "./precomputed_msa",
      "search_tool": "jackhmmer",
      "pairing_db": "uniprot",
      "pairing_db_fpath": "path_to_uniprot_database",
      "non_pairing_db_fpath": "path_to_mgnify_database",
      "msa_save_dir": "./searched_msa",
    },
  },
}
```
* `sequence`: A string representating a protein sequence, which can only contain the 20 standard amino acid type and X (UNK) for unknown residues.
* `count`: The number of copies of this protein chain (integer).
* `modifications`: An optional list of dictionaries that describe post-translational modifications.
  * `ptmType`: A string containing CCD code of the modification. 
  * `ptmPosition`: The position of the modified amino acid (integer).
* `msa`: A dictionary containing options for Multiple Sequence Alignment (MSA):
  * `precomputed_msa_dir`: The path to a directory containing precomputed MSAs. This directory should contain two specific files: "pairing.a3m" for MSAs used for pairing, and "non_pairing.a3m" for non-pairing MSAs. If you do not have precomputed MSAs, set this key to `None`.
  * `pairing_db`: The name of the genomic database used for pairing MSAs. Only "uniref100" and "uniprot" are supported. If you wish to precompute MSAs yourself, choose either the "uniref100" or "uniprot" database. However, if you lack precomputed MSAs and intend to search MSAs using our inference pipeline, only "uniprot" is available currently.
  * `search_tool`: The tool used for searching MSAs; currently, only "jackhmmer" is supported. This key is applicable only if precomputed MSAs are not available.
  * `pairing_db_fpath`: The path to the genomic database used for MSA pairing. This key is applicable only if precomputed MSAs are not available.
  * `non_pairing_db_fpath`: The path to the genomic database for non-pairing MSAs. This key is applicable only if precomputed MSAs are not available.
  * `msa_save_dir`: The directory to save the results of the MSA search. If not provided, searched MSAs will be removed after the inference process finishes. This key is applicable only if precomputed MSAs are not available.

##### dnaSequence
```json
{
  "dnaSequence": {
      "sequence": "GATTACA",
      "modifications": [
          {
              "modificationType": "CCD_6OG",
              "basePosition": 1
          },
          {
              "modificationType": "CCD_6MA",
              "basePosition": 2
          }
      ],
      "count": 1
  }
},
{
    "dnaSequence": {
        "sequence": "TGTAATC",
        "count": 1
    }
}
```
Please note that the `dnaSequence` type refers to a single stranded DNA sequence. If you
wish to model double-stranded DNA, please add a second `dnaSequence` entry representing
the sequence of the reverse complement strand.

* `sequence`: A string containing a DNA sequence; only letters A, T, G, C and N (unknown ribonucleotide) are allowed.
* `count`: The number of copies of this DNA chain (integer).
* `modifications`: An optional list of dictionaries describing of
the DNA chemical modifications:
  * `modificationType`: A string containing CCD code of modification.
  * `basePosition`: A position of the modified nucleotide (integer).

##### rnaSequence
```json
{
  "rnaSequence": {
      "sequence": "GUAC",
      "modifications": [
          {
              "modificationType": "CCD_2MG",
              "basePosition": 1
          },
          {
              "modificationType": "CCD_5MC",
              "basePosition": 4
          }
      ],
      "count": 1
  }
}
```
* `sequence`: A string representing the RNA sequence (single-stranded); only letters A, U, G, C and N (unknown nucleotides) are allowed.
* `count`: The number of copies of this RNA chain (integer).
* `modifications`: An optional list of dictionaries describing RNA chemical modifications:
  * `modificationType`: A string containing
    CCD code  of modification.
  * `basePosition`: The position of the modified nucleotide (integer).

##### ligand
```json
{
    "ligand": {
        "ligand": "CCD_ATP",
        "count": 1
    }
},
{
    "ligand": {
        "ligand": "FILE_your_file_path/atp.sdf",
        "count": 1
    }
},
{
    "ligand": {
        "ligand": "Nc1ncnc2c1ncn2[C@@H]1O[C@H](CO[P@@](=O)(O)O[P@](=O)(O)OP(=O)(O)O)[C@@H](O)[C@H]1O",
        "count": 1
    }
}
```
* `ligand`: A string representing the ligand. `ligand` can be one of the following three:
  * A string containing the CCD code of the ligand, prefixed with "CCD_". For glycans or similar structures, this can be a concatenation of multiple CCD codes, for example, "CCD_NAG_BMA_BGC".
  * A molecular SMILES string representing the ligand.
  * A path to a molecular structure file, prefixed with "FILE_", where the supported file formats are PDB, SDF, MOL, and MOL2. The file must include the 3D conformation of the molecule.

* `count` is the number of copies of this ligand (integer).

##### ion
```json
{
    "ion": {
        "ion": "MG",
        "count": 2
    }
},
{
    "ion": {
        "ion": "NA",
        "count": 3
    }
}
```
* `ion`: A string containing the CCD code for the ion. Note that, unlike ligands, the ion code **does not** start with "CCD_".
* `count`: The number of copies of this ion (integer).

#### covalent_bonds
```json
"covalent_bonds": [
            {
                "entity1": "2",
                "copy1": 1,
                "position1": "2",
                "atom1": "N6",
                "entity2": "3",
                "copy2": 1,
                "position2": "1",
                "atom2": "C1"
            }
]
```

The `covalent_bonds` section specifies covalent bonds between a polymer and a ligand, or between two ligands.
To define a covalent bond, two atoms involved in the bond must be identified. The following fields are used:

* `entity1`, `entity2`: The entity numbers for the two atoms involved in the bond. 
The entity number corresponds to the order in which the entity appears in the `sequences` list, starting from 1.
* `copy2`, `copy2`: The copy index (starting from 1) of the `left_entity` and `right_entity`, respectively. These fields are optional, but if specified, both `left_copy` and `right_copy` must be filled simultaneously or left empty at the same time. If neither field is provided, a bond will be created between all pairs of copies of the two entities. For example, if both entity1 and entity2 have two copies, a bond will be formed between entity1.copy1 and entity2.copy1, as well as between entity1.copy2 and entity2.copy2. In this case, the number of copies for both entities must be equal.
* `position1`, `position2` - The position of the residue (or ligand part) within the entity. 
The position value starts at 1 and can vary based on the type of entity:
  * For **polymers** (e.g., proteins, DNA, RNA), the position corresponds to the location of the residue in the sequence.
  * For **ligands** composed of multiple CCD codes, the position refers to the serial number of the CCD code.
  * For **single CCD code ligands**, or ligands defined by **SMILES** or **FILE**, the position is always set to 1.
    
* `atom1`, `atom2` - The atom names (or atom indices) of the atoms to be bonded.
  *  If the entity is a polymer or described by a CCD code, the atom names are consistent with those defined in the CCD.
  * If the entity is a ligand defined by SMILES or a FILE, atoms can be specified by their atom index. The atom index corresponds to the position of the atom in the file or in the SMILES string, starting from 0.

Deprecation Notice: The previous fields such as old `left_entity`, `right_entity`, and other fields starting with `left`/`right` have been updated to use `1` and `2` to denote the two atoms forming a bond. The current code still supports the old field names, but they may be deprecated in the future, leaving only the new field names.

### Format of the model output
The outputs will be saved in the directory provided via the `--dump_dir` flag in the inference script. The outputs include the predicted structures in CIF format and the confidence in JSON files. The `--dump_dir` will have the following structure:

```bash
├── <name>/  # specified in the input JSON file
│   ├── <seed>/  # specified via the `--seeds` flag in the inference script
│   │   ├── <name>_<seed>_sample_0.cif
│   │   ├── <name>_<seed>_summary_confidence_sample_0.json
│   │   └──... # the number of samples in each seed is specified via `--sample_diffusion.N_sample ` flag in the inference script
│   └──...
└── ...
```

The contents of each output file are as follows:
- `<name>_<seed>_sample_*.cif` - A CIF format text file containing the predicted structure
- `<name>_<seed>_summary_confidence_sample_*.json` - A JSON format text file containing various confidence scores for assessing the reliability of predictions. Here’s a description of each score:

    - `plddt` - Predicted Local Distance Difference Test (pLDDT) score. Higher values indicate greater confidence.
    - `gpde` - Globl Predicted Distance Error (PDE) score. Lower values indicate greater confidence.
    - `ptm` - Predicted TM-score (pTM). Values closer to 1 indicate greater confidence.
    - `iptm` - Interface Predicted TM-score, used to estimate the accuracy of interfaces between chains. Values closer to 1 indicate greater confidence.
    - `chain_ptm` - pTM score calculated for individual chains with the shape of [N_chains], indicating the reliability of specific chain structure.
    - `chain_pair_iptm`: Pairwise interface pTM scores between chain pairs with the shape of [N_chains, N_chains], indicating the reliability of specific chain-chain interactions.
    - `chain_iptm` - Average ipTM scores for each chain with the shape of [N_chains].
    - `chain_pair_iptm_global` - Averge `chain_iptm` between chain pairs with the shape of [N_chains, N_chains]. For interface containing a small molecule, ion, or bonded ligand chain (named `C*`), this value is equal to the `chain_iptm` value of `C*`. 
    - `chain_plddt` - pLDDT scores calculated for individual chains with the shape of [N_chains].
    - `chain_pair_plddt` - Pairwise pLDDT scores for chain pairs with the shape of [N_chains, N_chains].
    - `has_clash` - Boolean flag indicating if there are steric clashes in the predicted structure.
    - `disorder` - Predicted regions of intrinsic disorder within the protein, highlighting residues that may be flexible or unstructured.
    - `ranking_score` - Predicted confidence score for ranking complexes. Higher values indicate greater confidence.
    - `num_recycles`: Number of recycling steps used during inference.