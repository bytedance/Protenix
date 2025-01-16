## MSA data pipeline
If you download our released wwPDB dataset as in [training.md](./training.md), the mmcif_msa [450G] dir has the following directory structure.
  ```bash
  ├── seq_to_pdb_index.json [45M] # sequence to integers mapping file
  ├── mmcif_msa [450G] # msa files
    ├── 0
      ├── uniref100_hits.a3m
      ├── mmseqs_other_hits.a3m
    ├── 1
      ├── uniref100_hits.a3m
      ├── mmseqs_other_hits.a3m
    ├── 2
      ├── uniref100_hits.a3m
      ├── mmseqs_other_hits.a3m
    ...
    ├── 157201
      ├── uniref100_hits.a3m
      ├── mmseqs_other_hits.a3m

  ```

Each integer in the first-level directory under mmcif_msa (for example, 0, 1, 2, and 157201) represents a unique protein sequence. The key of `seq_to_pdb_index.json` is the unique protein sequence, and the value is the integer corresponding to the first-level subdirectory of mmcif_msa mentioned above.

This document is used to provide the steps to convert the MSA obtained from colabfold into the Protenix training format.

### Steps to get your own MSA data for training

#### Step1: get input protein sequence 
Run the following command:

```python
python3 scripts/msa/step1-get_prot_seq.py
```
you will get outputs in `scripts/msa/data/pdb_seqs` dir. The result dir is as follows,

```bash
  ├── pdb_index_to_seq.json # mapping integers to sequences
  ├── seq_to_pdb_index.json # mapping sequences to integers identifiers when saving MSA, This file is required in training for finding local MSA path from sequence
  ├── pdb_seq.fasta # Input of MSA
  ├── pdb_seq.csv # Intermediate Files
  ├── seq_to_pdb_id_entity_id.json # Intermediate Files
```

#### Step2: run msa search
We give detailed environment configuration and search commands in 

```python
scripts/msa/step2-get_msa.ipynb
```

The searched MSA is in `scripts/msa/data/mmcif_msa_initial`, The result dir is as follows,
```bash
  ├── 0.a3m
  ├── 1.a3m
  ├── 2.a3m
  ├── 3.a3m
  ├── pdb70_220313_db.m8
  ├── uniref_tax.m8 # record Taxonomy ID which is used by MSA Pairing
```
#### Steps3: MSA Post-Processing

The overall solution is to search the MSA containing taxonomy information only once for the unique sequence, and pair it according to the species information of each MSA. 

For MSA Post-Processing, Taxonomy ID from UniRef30 DB is added to MSAs and MSAs is split into `uniref100_hits.a3m` and `mmseqs_other_hits.a3m`, which correspond to `pairing.a3m` and `non_pairing.a3m` in inference stage respectively.

You can run:
```python
python3 scripts/msa/step3-uniref_add_taxid.py

python3 scripts/msa/step4-split_msa_to_uniref_and_others.py
```

The final pairing and non_pairing MSAs in `scripts/msa/data/mmcif_msa` is as follows:


```
>query
GPTHRFVQKVEEMVQNHMTYSLQDVGGDANWQLVVEEGEMKVYRREVEENGIVLDPLKATHAVKGVTGHEVCNYFWNVDVRNDWETTIENFHVVETLADNAIIIYQTHKRVWPASQRDVLYLSVIRKIPALTENDPETWIVCNFSVDHDSAPLNNRCVRAKINVAMICQTLVSPPEGNQEISRDNILCKITYVANVNPGGWAPASVLRAVAKREYPKFLKRFTSYVQEKTAGKPILF
>UniRef100_A0A0S7JZT1_188132/	246	0.897	6.614E-70	2	236	237	97	331	332
--THRFADKVEEMVQNHMTYSLQDVGGDANWQLVIEEGEMKVYRREVEENGIVLDPLKATHAVKGVTGHEVCHYFWDTDVRNDWETTIDNFNVVETLSDNAIIVYQTHKRVWPASQRDILFLSAIRKILAKNENDPDTWLVCNFSVDHDKAPPTNRCVRAKINVAMICQTLVSPPEGDKEISRDNILCKITYVANVNPGGWAPASVLRAVAKREYPKFLKRFTSYVQEKTAGNPILF
>UniRef100_A0A4W6GBN4_8187/	246	0.893	9.059E-70	2	236	237	373	607	608
--THRFANKVEEMVQNHMTYSLQDVGGDANWQLVIEEGEMKVYRREVEENGIVLDPLKATHSVKGVTGHEVCHYFWDTDVRMDWETTIENFNVVEKLSENAIIVYQTHKRVWPASQRDVLYLSAIRKIMATNENDPDTWLVCNFSVDHNNAPPTNRCVRAKINVAMICQTLVSPPEGDKEISRDNILCKITYVANVNPGGWAPASVLRAVAKREYPKFLKRFTSYVQEKTAGKPILF
```

```
>query
MAEVIRSSAFWRSFPIFEEFDSETLCELSGIASYRKWSAGTVIFQRGDQGDYMIVVVSGRIKLSLFTPQGRELMLRQHEAGALFGEMALLDGQPRSADATAVTAAEGYVIGKKDFLALITQRPKTAEAVIRFLCAQLRDTTDRLETIALYDLNARVARFFLATLRQIHGSEMPQSANLRLTLSQTDIASILGASRPKVNRAILSLEESGAIKRADGIICCNVGRLLSIADPEEDLEHHHHHHHH
>MGYP001165762451	218	0.325	1.019E-59	5	230	244	3	228	230
-----DKVEFLKGVPLFSELPEAHLQSLGELLIERSYRRGATIFFEGDPGDALYIVRSGIVKISRVAEDGREKTLAFLGKGEPFGEMALIDGGPRSAIAQALEATSLYALHRADFLAALTENPALSLGVIKVLSARLQQANAQLMDLVFRDVRGRVAQALLDLARR-HGVPLTNGRMISVKLTHQEIANLVGTARETVSRTFAELQDSGIIRIeGRNIVLLDAAQLEGYAAG-------------
>A0A160T8V6	218	0.285	1.019E-59	0	227	244	0	229	237
MPTTRDsnAVQALQVVPFFANLPEDHVAALAKALVPRRFSPGQVIFHLGDPGGLLYLISRGKIKISHTTSDGQEVVLAILGPGDFFGEMALIDDAPRSATAITLEPSETWTLHREEFIQYLTDNPEFALHVLKTLARHIRRLNTQLADIFFLDLPGRLARTLLNLADQ-YGRRAADGTIIDLSLTQTDLAEMTGATRVSINKALGRFRRAGWIQvTGRQVTVLDRAALEAL----------------
>AP58_3_1055460.scaffolds.fasta_scaffold1119545_2	216	0.304	3.581E-59	10	225	244	5	221	226
----------LSRVPLFAELPPERIHELAQSVRRRTYHRGETIFHKGDPGNGLYIIAAGQVKIVLPSEMGEEAMLAVLEGGEFFGELALFDGLPRSATVVAVQNAEVLVLHRDDFMSFVGRNPEVVSALFAALSRRLRDADEMIEDAIFLDVPGRLAKRLLDLAEKHGRAEEKGGVAIDLKLTQQDLAAMVGATRESVNKHLGWMRDHGLIQLDRqRIVILKPDDLR------------------
```
### Format of MSA 
In `uniref100_hits.a3m`(training stage) or `pairing.a3m`(inference stage), the header must starts with the following format, which we use for pairing:
```
>UniRef100_{hitname}_{taxonomyid}/
```

we also provide a pipeline of local Colabfold_search to Generate Protenix-Compatible MSAs in [colabfold_compatiable_msa.md](./colabfold_compatiable_msa.md).