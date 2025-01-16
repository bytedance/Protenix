### Using Local Colabfold_search to Generate Protenix-Compatible MSA

Colabfold provides an easy-to-use and efficient MSA search pipeline that's ideal for generating MSAs during inference. Unfortunately, this pipeline cannot fully match Protenix's MSA search process designed for training, as the current `colabfold_search` omits species information in the MSA, preventing correct pairing by Protenix's data pipeline. To address this issue, we provide the `scripts/colabfold_msa.py` script, which post-processes `colabfold_search` results by adding pseudo taxonomy IDs to paired MSAs to match Protenix's data pipeline.

Here's an example:
```bash
python3 scripts/colabfold_msa.py examples/dimer.fasta <path/to/colabfold_db> dimer_colabfold_msa --db1 uniref30_2103_db --db3 colabfold_envdb_202108_db --mmseqs_path <path/to/mmseqs> 
```

#### Configuring Colabfold_search
Installation of colabfold and mmseqs2 is required.

colabfold can be installed with: `pip install colabfold[alphafold]`. 

Build MMseqs2 from source:

```bash
wget https://github.com/soedinglab/MMseqs2/archive/refs/tags/16-747c6.tar.gz
tar xzf 16-747c6.tar.gz 
cd MMseqs2-16-747c6/
mkdir build && cd build
cmake -DCMAKE_BUILD_TYPE=RELEASE -DCMAKE_INSTALL_PREFIX=. ..
make -j8
make install
```

Download ColabFold database:
```bash
git clone https://github.com/sokrypton/ColabFold.git
cd ColabFold
# Configure database:
MMSEQS_NO_INDEX=1 ./setup_databases.sh <path/to/colabfold_db>
```