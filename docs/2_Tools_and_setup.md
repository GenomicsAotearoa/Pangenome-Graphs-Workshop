# 2. Tools & setup

!!! screwdriver-wrench "Tools used for the pangenome graph pipeline"

    - Graph construction using the PanGenome Graph Builder (PGGB) (https://github.com/pangenome/pggb)
    
    - Graph manipulation through the Optimized Dynamic Genome/Graph Implementation (ODGI)(https://github.com/pangenome/odgi)
    
    - Variant calling for Next-Generation Sequencing (NGS) data using the VG toolkit(https://github.com/vgteam/vg)

    - Circlator (https://sanger-pathogens.github.io/circlator/)
    - Mash (https://github.com/marbl/Mash)
    - SAMtools (https://github.com/samtools/samtools)
    - BCFtools (https://github.com/samtools/bcftools)
    - simuG (https://github.com/yjx1217/simuG)
    - PGGE (https://github.com/pangenome/pgge)

!!! info "Running the pggb pipeline on NeSI"

    PGGB, ODGi, VG, circulator, Mash, SAMTools et.al. have been installed on NeSI as modules. We need to load each module first.  
    ```bash
    module load pggb
    module load SAMtools
    module load Circlator/1.5.5-gimkl-2022a-Python-3.10.5
    module load Mash/2.3-GCC-11.3.0
    ......
    ```

??? info "Running the pggb workflow locally"

    ### Obtaining pggb
    From https://github.com/pangenome/pggb, you can find the details about installing pggb with Docker, Singularity, bioconda, guix, or by manually building its dependencies.
    
    ### Using pggb via Singularity
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
    Finally, run pggb from the Singularity image. For Singularity to be able to read and write files to a directory on the host     operating system, we need to 'bind' that directory using the -B option and pass the pggb command as an argument.
    ```bash
    singularity run -B ${PWD}/data:/data ../pggb_latest.sif pggb -i /data/HLA/DRB1-3123.fa.gz -p 70 -s 3000 -G 2000 -n 10 -t 16 -v -V 'gi|568815561:#' -o /data/out -M -m
    ```

!!! info ""

    ### *Neisseria meningitidis* data set
    
    For this workshop, we utilized the genomes of the bacterium _Neisseria meningitidis_ as a representative example.
    _Neisseria (N.) meningitidis_, also known as the meningococcus pathogen, is the primary agent responsible for invasive meningococcal diseases such as meningitis and septicemia, causing isolated incidents, outbreaks, and epidemics worldwide. The genome of this bacterium spans approximately 2.1 to 2.4 Mb and possesses a GC content ranging from 51-52%. One striking characteristic of _N. meningitidis_ genomes is their high recombination rate, which largely fuels the extensive genetic diversity     within this species. In this workshop, we utilized five genome assemblies of _N. meningitidis_ to assess the pangenome pipeline, covering pangenome graph construction to variant calling. 

| genomes                             | ASM IDs   |GCF IDs    | SEROGROUP  | Sequence type | Clonal Complex   |
|:-----                               |----------:|----------:|-----------:|--------------:|-----------------:|
|NC_017518.1 Neisseria meningitidis NZ-05/33	| ASM19152v1	|GCF_000191525.1	|B	|42	|ST-41/44	|
|NC_003112.2 Neisseria meningitidis MC58    	| ASM880v1	|GCF_000008805.1 	|-	|74	|ST-32	|
|NZ_CP007668.1 Neisseria meningitidis M0579  	| ASM102983v1|GCF_001029835.1		|B	|-	|ST-41/44	|
|NZ_CP016880.1 Neisseria meningitidis strain M07165	| ASM170367v1	|GCF_001703675.1	|W	|11	|ST-11	|
|NZ_CP020423.2 Neisseria meningitidis strain FDAARGOS_212	| ASM207367v2	|GCF_002073675.2	|C	|-	|ST16521	|

??? info "How the *Neisseria meningitidis* genomes were formatted for this workshop"

    The following code **DOES NOT** need to be run, but is provided here to show how the *Neisseria meningitidis* genomes
    were downloaded and prepared for analysis.
    
    ```bash
    # Create a new directory called nm_genomes and change to that directory
    mkdir nm_genomes
    cd nm_genomes
    ```
    
    Download the genome assemblies from NCBI and uncompress.  In the Unix environment you can use the `curl` command.

    ```bash
    #NC_017518.1
    curl -OJX GET "https://api.ncbi.nlm.nih.gov/datasets/v2alpha/genome/accession/GCF_000191525.1/download?include_annotation_type=GENOME_FASTA,GENOME_GFF,RNA_FASTA,CDS_FASTA,PROT_FASTA,SEQUENCE_REPORT&filename=GCF_000191525.1.zip" -H "Accept: application/zip"

    # ACC NUM?
    curl -OJX GET "https://api.ncbi.nlm.nih.gov/datasets/v2alpha/genome/accession/GCF_001029835.1/download?include_annotation_type=GENOME_FASTA,GENOME_GFF,RNA_FASTA,CDS_FASTA,PROT_FASTA,SEQUENCE_REPORT&filename=GCF_001029835.1.zip" -H "Accept: application/zip"

    # ACC NUM?
    curl -OJX GET "https://api.ncbi.nlm.nih.gov/datasets/v2alpha/genome/accession/GCF_001703675.1/download?include_annotation_type=GENOME_FASTA,GENOME_GFF,RNA_FASTA,CDS_FASTA,PROT_FASTA,SEQUENCE_REPORT&filename=GCF_001703675.1.zip" -H "Accept: application/zip"

    # ACC NUM?
    curl -OJX GET "https://api.ncbi.nlm.nih.gov/datasets/v2alpha/genome/accession/GCF_002073675.2/download?include_annotation_type=GENOME_FASTA,GENOME_GFF,RNA_FASTA,CDS_FASTA,PROT_FASTA,SEQUENCE_REPORT&filename=GCF_002073675.2.zip" -H "Accept: application/zip"

    # ACC NUM?
    curl -OJX GET "https://api.ncbi.nlm.nih.gov/datasets/v2alpha/genome/accession/GCF_000008805.1/download?include_annotation_type=GENOME_FASTA,GENOME_GFF,RNA_FASTA,CDS_FASTA,PROT_FASTA,SEQUENCE_REPORT&filename=GCF_000008805.1.zip" -H "Accept: application/zip"
    ```

    On NeSI, a slurm job can be run to process the .fna genomes.  The contents of the slurm job file (`unzip_genomes.sl`) 
    are as follows: 

    ```bash
    #!/usr/bin/bash

    #SBATCH --account       nesi02659
    #SBATCH --job-name      extract_fna
    #SBATCH --cpus-per-task 8
    #SBATCH --mem           4G
    #SBATCH --time          1:00:00

    data=$HOME/nm_genomes/*.zip

    for f in $data

    do

        x=$(basename $f .zip)
        echo ${x}

        unzip $x.zip

        cp $HOME/nm_genomes/ncbi_dataset/data/${x}/*_genomic.fna /$HOME/nm_genomes/

        rm -rf ncbi_dataset

    done
    ```

    The slurm job can be run via:

    ```
    sbatch unzip_genomes.sl
    ```

    Remove unneeded files:
    
    ```bash
    rm cds_from_genomic.fna
    rm *.zip
    rm README.md
    rm slurm*.out
    ```

    Use the `cat` command to combine genomes into a single fasta file:

    ```
    cat *_genomic.fna > 5NM.fa
    ```

### Setting up your project directory and download the datasets

!!! terminal "code"

    ```bash
    # Create a new directory in somewhere and change to that directory
    mkdir pg_workshop
    cd pg_workshop
    # Keep a note of the absolute path of your directory
    pwd
    ```
    !!! success "Output"

        ```
        /home/<YOUR_USER_ID>/pg_workshop
        ```
    
    ```bash
    # Downloading and preparing datasets
    git clone https://github.com/ZoeYang2020/dataset_for_pg_workshop
    
    # Copy the 5NM.fa dataset to your work directory
    cp /home/<YOUR_USER_ID>/datasets_for_PangenomeGraphConstruction_pg_workshop/5NM.fa ./
    ```
