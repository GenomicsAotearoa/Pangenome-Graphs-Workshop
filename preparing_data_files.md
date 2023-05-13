# Preparing Required Data Files for the exersice with simuG
[simuG](https://github.com/yjx1217/simuG) is a light-weighted tool which can be used to simulate sequences with pre-defined number of random genomic variants against a reference sequence (Please refer the [Git Page](https://github.com/yjx1217/simuG) for more details). 

### 1. Setting up the tool 
Please refer the installation instruction in the [Git Page](https://github.com/yjx1217/simuG).
```bash
git clone https://github.com/yjx1217/simuG.git
cd simuG
perl simuG.pl -h
perl vcf2model.pl -h
```

### 2. Data Source
You can download the reference sequence (Neisseria gonorrhoeae FA 1090 genome assembly) from the site [Genome assembly ASM684v1](https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_000006845.1/) and get the file \ncbi_dataset\data\GCA_000006845.1\GCA_000006845.1_ASM684v1_genomic.fna from he zip file ncbi_dataset.zip. Rename it as "ref.fa" to refere easily in next steps. 

