## MSA Data Pipeline

This document describes the MSA (Multiple Sequence Alignment) data organization in Protenix and provides a step-by-step guide to generating your own compatible MSA data.

### Overview of MSA Data Structure

If you are using the released wwPDB dataset as described in [training_inference_instructions.md](./training_inference_instructions.md), your data directory should be organized as follows:

```bash
├── common
│   ├── seq_to_pdb_index.json # Mapping from sequence to directory ID
│   └── ...                   # Other metadata (CCD, release dates, etc.)
├── mmcif_msa_template                 # Root directory for MSA files
│   ├── 0
│   │   ├── pairing.a3m       # MSAs paired by taxonomy ID
│   │   ├── non_pairing.a3m   # Unpaired MSAs
│   │   └── hmmsearch.a3m     # Templates from HMMER search
│   ├── 1
│   │   ├── pairing.a3m
│   │   ├── non_pairing.a3m
│   │   └── hmmsearch.a3m
│   └── ...

```

Each unique protein sequence in the dataset is assigned an integer index (e.g., `0`, `1`, `2`). These indices correspond to the subdirectories under `mmcif_msa_template`. The `seq_to_pdb_index.json` file is used by the data loader to resolve a sequence to its corresponding directory ID.

---

### Generating Custom MSA Data for Training

If you are working with a new dataset, follow these steps to generate compatible MSA data.

#### Prerequisites

- Ensure all input mmCIF files are stored in a single directory (e.g., `scripts/msa/data/mmcif`).
- Python environment with required dependencies installed.

#### Step 1: Extract Unique Protein Sequences

Run the following script to extract sequences from mmCIF files and assign index IDs:

```bash
python3 scripts/msa/step1-get_prot_seq.py --cif_dir <path_to_cif_dir> --out_dir scripts/msa/data/pdb_seqs
```

**Key Outputs in `scripts/msa/data/pdb_seqs`:**
- `pdb_seq.fasta`: Combined FASTA file containing all unique protein sequences, used as input for the MSA search.
- `seq_to_pdb_index.json`: Mapping from sequence to index ID. **Move this file to your dataset root after generation.**
- `seq_to_pdb_id_entity_id.json`: Mapping from sequence to PDB/entity IDs, used in Step 3.
- `pdb_index_to_seq.json`: Reverse mapping from index ID to sequence.

#### Step 2: Execute MSA Search

Protenix requires MSAs with taxonomic information for pairing. You can use the provided notebook as a template for the search process:

```python
scripts/msa/step2-get_msa.ipynb
```

The raw search results should be saved in `scripts/msa/data/mmcif_msa_initial` (default format is `.a3m` files named by index ID, e.g., `0.a3m`).

#### Step 3: MSA Post-Processing

This step converts raw A3M files into the final Protenix-compatible format by appending taxonomy IDs and splitting hits.

1. **Append Taxonomy IDs**:
   ```bash
   python3 scripts/msa/step3-uniref_add_taxid.py
   ```
2. **Split and Organize**:
   ```bash
   python3 scripts/msa/step4-split_msa_to_uniref_and_others.py --input_msa_dir scripts/msa/data/mmcif_msa_with_taxid --output_msa_dir scripts/msa/data/mmcif_msa
   ```

The processed MSAs will be saved in `scripts/msa/data/mmcif_msa`, organized by their sequence index.

---

### MSA Header Format Requirements

To enable correct taxonomic pairing, the headers in `pairing.a3m` must follow specific patterns. Protenix supports both UniRef and UniProt header formats.

**UniRef Format:**
```text
>UniRef100_{hit_name}_{taxonomy_id}/
```

**UniProt Format:**
```text
>tr|ACCESSION|ID_SPECIES/START-END [optional description] OS=Species Name OX=TaxonomyID ...
```

The system uses the following regex patterns to extract `SpeciesId`:

```python
_UNIPROT_REGEX = re.compile(
    r"(?:tr|sp)\|[A-Z0-9]{6,10}(?:_\d+)?\|(?:[A-Z0-9]{1,10}_)(?P<SpeciesId>[A-Z0-9]{1,5})"
)
_UNIREF_REGEX = re.compile(r"^UniRef100_[^_]+_([^_/]+)")
```

### Next Steps

For a simplified local search pipeline using Colabfold, please refer to [colabfold_compatible_msa.md](./colabfold_compatible_msa.md).

---

## Template Search

Template search is used to find structural templates for protein sequences. It uses `hmmsearch` to search against the `pdb_seqres` database based on the profile of existing MSAs.

### Usage

The template search can be executed via:

```bash
python3 runner/template_search.py
```

The core function `run_template_search` performs the following:
- Takes a directory containing `pairing.a3m` and `non_pairing.a3m` as input.
- Builds an HMM profile from these MSAs and searches the database.
- Outputs a merged `hmmsearch.a3m` file in the same directory.

### Database

The default database is `pdb_seqres_2022_09_28.fasta`. The script will automatically download it from ByteDance TOS if it is missing from `PROTENIX_ROOT_DIR/search_database`.

---

## RNA MSA Search

Protenix provides a specialized pipeline for RNA MSA search using `nhmmer`.

### Usage

To run RNA MSA search for a specific sequence:

```bash
python3 runner/rna_msa_search.py
```

The pipeline searches across three major RNA databases. You can use the clustered versions of the databases released by AlphaFold 3 or provide your own custom databases:
1. **Rfam**: `rfam_14_9_clust_seq_id_90_cov_80_rep_seq.fasta`
2. **RNAcentral**: `rnacentral_active_seq_id_90_cov_80_linclust.fasta`
3. **NT-RNA**: `nt_rna_2023_02_23_clust_seq_id_90_cov_80_rep_seq.fasta`

### Output

The search results from all databases are merged, deduplicated, and saved as `rna_msa.a3m`. This file should be referenced in the inference JSON as `unpairedMsaPath` for RNA chains.

### Prerequisites

- **HMMER**: Ensure `hmmer` is installed and binaries (`nhmmer`, `hmmalign`, `hmmbuild`) are in your PATH.
- **Storage**: RNA databases are large; ensure several gigabytes of free space in your search database directory.


