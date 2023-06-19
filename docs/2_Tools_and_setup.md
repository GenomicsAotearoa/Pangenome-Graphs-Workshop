# Tool & setup

## Tools used in this pipeline


!!! screwdriver-wrench "Tools used for the pangenome graph pipeline"

    - Graph construction using the PanGenome Graph Builder (PGGB) (https://github.com/pangenome/pggb)
    
    - Graph manipulation through the Optimized Dynamic Genome/Graph Implementation (ODGI)(https://github.com/pangenome/odgi)
    
    - Variant calling for Next-Generation Sequencing (NGS) data using the VG toolkit(https://github.com/vgteam/vg)

    - circlator (https://sanger-pathogens.github.io/circlator/)
    - Mash (https://github.com/marbl/Mash)
    - SAMtools (https://github.com/samtools/samtools)
    - BCFtools (https://github.com/samtools/bcftools)
    - simuG (https://github.com/yjx1217/simuG)
    - PGGE (https://github.com/pangenome/pgge)

## Running on NeSi
PGGB, ODGi, VG, circulator, Mash, SAMTools et.al. have been installed on Nesi as modules. We need to load each module first. 
```bash
module load pggb
module load SAMtools
module load Circlator/1.5.5-gimkl-2022a-Python-3.10.5
module load Mash/2.3-GCC-11.3.0
......

```
## Running locally
From https://github.com/pangenome/pggb, you can find the details about installing pggb with Docker, Singularity, bioconda, guix, or by manually building its dependencies.
### Singularity
Many managed HPCs utilize Singularity as a secure alternative to docker. Fortunately, docker images can be run through Singularity seamlessly.
First pull the docker file and create a Singularity SIF image from the dockerfile. This might take a few minutes.
```bash
singularity pull docker://ghcr.io/pangenome/pggb:latest
```
Next clone the pggb repo and cd into it
```bash
git clone --recursive https://github.com/pangenome/pggb.git
cd pggb
```
Finally, run pggb from the Singularity image. For Singularity to be able to read and write files to a directory on the host operating system, we need to 'bind' that directory using the -B option and pass the pggb command as an argument.
```bash
singularity run -B ${PWD}/data:/data ../pggb_latest.sif pggb -i /data/HLA/DRB1-3123.fa.gz -p 70 -s 3000 -G 2000 -n 10 -t 16 -v -V 'gi|568815561:#' -o /data/out -M -m
```



# Data
For this workshop, we utilized the genomes of the bacterium Neisseria meningitidis as a representative example.
Neisseria (N.) meningitidis, also known as the meningococcus pathogen, is the primary agent responsible for invasive meningococcal diseases such as meningitis and septicemia, causing isolated incidents, outbreaks, and epidemics worldwide. The genome of this bacterium spans approximately 2.1 to 2.4 Mb and possesses a GC content ranging from 51-52%. One striking characteristic of N. meningitidis genomes is their high recombination rate, which largely fuels the extensive genetic diversity within this species. In this workshop, we utilized genome assemblies of N. meningitidis to assess the pangenome pipeline, covering pangenome graph construction to variant calling. 


## Setting up your project directory


# Create a new directory in somewhere and change to that directory
mkdir pg_workshop
cd pg_workshop
# Keep a note of the absolute path of your directory
pwd
/home/$youraccount/pg_workshop
```

## Preparing input data

### Genome Availability 
The National Library of Medicine is the largest library focused on biomedicine worldwide, serving as the central hub for biomedical informatics and computational biology. It has many genome assembly data and [Genome assembly ASM19152v1](https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_000191525.1/) will be used for this workshop. 

---
# Procedure 
### 1. Downloading and preparing assembly data file 4Sim.fa

Please follow the procedure described in this [page](./preparing_data_files.md)

### 2. Creating an index for the sequence file and check

!!! terminal "code"

    ```bash
    # Use SAMtools to create the index file
    # In NeSI environment you will have to load the command first
    
    module load SAMtools
    
    samtools faidx ASM19152v1_pgsim.fa 
    
    cat ASM19152v1_pgsim.fa.fai
    ```

!!! success "Output"

    ```
    NC_017518.1     2248966 64      80      81
    NC_017518.1_SNP_5000    2248966 2277165 2248966 2248967
    NC_017518.1_INDEL_5000  2249048 4526156 2249048 2249049
    NC_017518.1_SNP_4000_INDEL_4000 2242147 6775238 2242147 2242148
    NC_017518.1_SNP_4000_INDEL_4000_INV_4   2242147 9017425 2242147 2242148
    NC_017518.1_SNP_4000_INDEL_4000_CNV_4   2415498 11259612        2415498 2415499
    ```
    
As per the index this assembly consists of 6 samples described in the below table. 

| Name                                | Length    | SNPs   | INDELs | INV | CNV |
|:-----                               |----------:|-------:|-------:|----:|----:|
|NC_017518.1 (Reference)              | 2,248,966 |     N/A|     N/A| N/A |  N/A|
|NC_017518.1_SNP_5000                 | 2,248,966 |   5,000|       0|   0 |   0 |
|NC_017518.1_INDEL_5000               | 2,249,048 |       0|   5,000|   0 |   0 |
|NC_017518.1_SNP_4000_INDEL_4000      | 2,153,883 |   4,000|   4,000|   0 |   0 |
|NC_017518.1_SNP_4000_INDEL_4000_INV_4| 2,242,147 |   4,000|   4,000|   4 |   0 |
|NC_017518.1_SNP_4000_INDEL_4000_CNV_4| 2,415,498 |   4,000|   4,000|   0 |   4 |
