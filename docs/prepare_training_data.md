# Start with CIF files and prepare your own training data.

## Data Preparation

1. **Prepare CIF Files**: Place the CIF files you want to convert into training data in a folder. Alternatively, you can use a `txt` file to record the paths to these CIF files, with each line corresponding to the path of a specific CIF file.

2. **Prepare Protein Clustering File (Optional)**: The protein clustering file contains category information for each `[PDB ID]_[Entity ID]`. In the Protenix training data, we cluster protein sequences using a 40% sequence identity threshold.

   You can download the official clustering results file provided by RCSB PDB using the following command, and use it directly:
   ```bash
   wget https://cdn.rcsb.org/resources/sequence/clusters/clusters-by-entity-40.txt
   ```

    If you prefer to perform your own clustering of protein sequences, ensure the final results are formatted as a text file like this:
    Each line represents a cluster, containing `[PDB ID]_[Entity ID]` entries separated by spaces.

3. **Update the CCD (Chemical Component Dictionary) Cache File (If needed)**: We provide a pre-processed file, with a cutoff date of 2024-06-08, that records the reference conformers for each CCD Code. If the training data you're preparing is more recent than this date, there may be issues with some CCD Codes might be missing. For example, the CCD Code "WTB," appearing in the PDB ID: 8P3K released on 2024-11-20, is not defined in the previously provided CCD file. In such cases, you need to run the following script to download and update the CCD CIF files:

    ```bash
    python3 scripts/gen_ccd_cache.py -c [ccd_cache_dir] -n [num_cpu]
    ```

    After running the script, three files will be generated in the specified "ccd_cache_dir":
    
    - `components.cif` (CCD CIF file downloaded from RCSB)
    - `components.cif.rdkit_mol.pkl` (pre-processed dictionary, where the key is the CCD Code and the value is an RDKit Mol object with 3D structure)
    - `components.txt` (a list containing all the CCD Codes)

    When running Protenix, it first uses 
    ```bash
    `release_data/ccd_cache/components.cif`
    `release_data/ccd_cache/components.cif.rdkit_mol.pkl`
    ```
    if unavailable, it switches to 
    ```bash
    `release_data/ccd_cache/components.v20240608.cif`
    `release_data/ccd_cache/components.v20240608.cif.rdkit_mol.pkl`
    ```
    Notes:
    - The `-c` parameter is optional. If not specified, files will be saved in the "release_data/ccd_cache" folder within the Protenix code directory by default.
    - You can add the `-d` parameter when running the script to skip the CIF file download step, in which case the script will directly process the "components.cif" file located in the "ccd_cache_dir".

## Data Preprocessing
Execute the script to preprocess the data:
```bash
python3 scripts/prepare_training_data.py -i [input_path] -o [output_csv] -b [output_dir] -c [cluster_txt] -n [num_cpu]
```

The preprocessed structures will be saved as `.pkl.gz` files. Additionally, a `CSV` file will be generated to catalog the chains and interfaces within these structures, which will facilitate sampling during the training process.

You can view the explanation of the parameters by using the `--help` command.
```
python3 scripts/prepare_training_data.py --help
```

Note that there is an optional parameter `-d` in the script. When this parameter is not used, the script processes CIF files downloaded from RCSB PDB by applying the full set of WeightedPDB training data filters. These filters include:

- Removing water molecules
- Removing hydrogen atoms
- Deleting polymer chains composed entirely of unknown residues
- Eliminating chains where the Cα distance between adjacent numbered residues exceeds 10 angstroms
- Removing elements labeled as "X"
- Deleting chains where no residues have been resolved
- When the number of chains exceeds 20, selecting one central atom from those capable of forming interfaces and retaining the 20 nearest chains to it. If a ligand is covalently bonded to a polymer, it is considered as one chain together. Additionally, if the number of chains is greater than 20 but the total number of tokens in these chains is less than 5120, more chains will be retained until the 5120 token limit is reached.
- Removing chains with one-third of their heavy atoms colliding

For CIF files generated through model inference where these filtering steps aren't desired, you can run the script with the `-d` parameter, which disables all these filters. The CIF structure will not be expanded to Assembly 1 in this case.


## Output Format
### Bioassembly Dict
In the folder specified by the `-b` parameter of the data preprocessing script, a corresponding `[pdb_id].pkl.gz` file is generated for each successfully processed CIF file. This file contains a dictionary saved with `pickle.dump`, with the following contents:
```
| Key                        | Value Type    | Description                                                                   |
|----------------------------|---------------|-------------------------------------------------------------------------------|
| pdb_id                     | str           | PDB Code                                                                      |
| assembly_id                | str           | Assembly ID                                                                   |
| sequences                  | dict[str, str]| Key is polymer's label_entity_id, value is canonical_sequence                 |
| release_date               | str           | PDB's Release Date                                                            |
| num_assembly_polymer_chains| int           | Number of assembly polymer chains (pdbx_struct_assembly.oligomeric_count)     |
| num_prot_chains            | int           | Number of protein chains in AtomArray                                         |
| entity_poly_type           | dict[str, str]| Key is polymer's label_entity_id, value is corresponding to entity_poly.type  |
| resolution                 | float         | Resolution; if no resolution, value is -1                                     |
| num_tokens                 | int           | Number of tokens                                                              |
| atom_array                 | AtomArray     | AtomArray from structure processing                                           |
| token_array                | TokenArray    | TokenArray generated based on AtomArray                                       |
| msa_features               | None          | (Placeholder)                                                                 |
| template_features          | None          | (Placeholder)                                                                 |
```

### Indices CSV
After the script successfully completes, a CSV file will be generated in the directory specified by `-o`. 
Each row contains information about a pre-processed chain or interface, and the content of each column is described as follows: 
```
| Column Name    | Value Type | Meaning                                                                | Required |
|----------------|------------|------------------------------------------------------------------------|----------|
| type           | str        | "chain" or "interface"                                                 | Y        |
| pdb_id         | str        | PDB Code (entry.id)                                                    | Y        |
| cluster_id     | str        | Cluster_id of the chain/interface                                      | Y        |
| assembly_id    | str        | Assembly id                                                            | N        |
| release_date   | str        | Release date                                                           | N        |
| resolution     | float      | Resolution; if no resolution, value is -1                              | N        |
| num_tokens     | int        | Number of tokens in AtomArray of Bioassembly Dict                      | N        |
| num_prot_chains| int        | Number of protein chains in AtomArray of Bioassembly Dict              | N        |
| eval_type      | str        | Classification used for evaluation                                     | N        |
| entity_1_id    | str        | Chain 1's label_entity_id                                              | Y        |
| chain_1_id     | str        | Chain 1's chain ID                                                     | Y        |
| mol_1_type     | str        | Chain 1's corresponding mol_type ("protein", "nuc", "ligand", "ions")  | Y        |
| sub_mol_1_type | str        | Sub-classification of Chain 1's entity corresponding to mol_type       | N        |
| cluster_1_id   | str        | Chain 1's cluster ID                                                   | Y        |
| entity_2_id    | str        | Chain 2's label_entity_id                                              | Y        |
| chain_2_id     | str        | Chain 2's chain ID                                                     | Y        |
| mol_2_type     | str        | Chain 2's corresponding mol_type ("protein", "nuc", "ligand", "ions")  | Y        |
| sub_mol_2_type | str        | Sub-classification of Chain 2's entity corresponding to mol_type       | N        |
| cluster_2_id   | str        | Chain 2's cluster_id                                                   | Y        |
```
Notes: 
- In the table, columns marked with 'Y' under 'Required' indicate that these columns are essential for training. If you are creating your own CSV for training purposes, these columns must be included. Columns marked with 'N' are optional and can be excluded. 
- For rows where the "type" is "chain", the values in columns related to Chain 2 should all be filled with empty strings.
