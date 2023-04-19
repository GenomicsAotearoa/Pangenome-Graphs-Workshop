# Pangenome Graphs Workshop
<p align="justify">
This repository contains material to construct pangenome graphs for a small bacteria dataset detailing every step in the Nesi environment. This study includes an anaysis of Neisseria Bacteria genome sequence data with 4Sim data samples to construct pangenome graphs to identify genetic variation and structural variance.
</p>

### Brief Overview of Pangenome Graphs
<p align="justify">
A pangenome graph is a graphical representation of the collective genomic information of a set of related organisms. Unlike a traditional linear genome assembly, which represents a single consensus genome sequence, a pangenome graph captures the genetic variation and structural diversity within a population. Pangenome graphs are constructed by integrating multiple genome sequences into a single graph structure that represents the entire set of genetic elements present in the population. This graph structure allows for the identification and visualization of genomic variations, such as single nucleotide polymorphisms (SNPs), insertions, deletions, and structural variations, as well as the relationships between different genomic regions. Pangenome graphs have become an important tool in genomics research, especially in the study of bacterial and viral populations, where genetic diversity is often high.
 </p>

### Neisseria Genome Assembly
<p align="justify">
Neisseria is a genus of Gram-negative bacteria that are typically found in the mucous membranes and epithelial tissues of humans and animals. There are various species of Neisseria bacteria, some of which are harmless commensals, while others can cause diseases. Two well-known pathogenic species are Neisseria meningitidis, which can cause meningitis and septicemia, and Neisseria gonorrhoeae, which is the causative agent of the sexually transmitted infection gonorrhea. Other species of Neisseria are generally considered harmless and reside as commensals in the oral and/or nasopharynx of humans.
 </p>

### PanGenome Graph Builder [(PGGB)](https://github.com/pangenome/pggb)
<p align="justify">
To investigate and analyse Neisseria Bacteria the PanGenome Graph Builder has been used. PanGenome Graph Builder (PGGB) is a computational tool used in genomics research to construct pan-genome graphs from large sets of genomic data. A pan-genome is the collection of all the genes and non-coding sequences present in a species or a group of related organisms. PGGB constructs a pan-genome graph, which is a data structure that represents the entire set of genes and genetic variations in a population. This graph can be used to study genetic diversity, gene function, and evolution. PGGB is designed to be scalable and efficient, making it suitable for large-scale genomic analyses. It is an open-source tool that can be used freely by researchers in the field of genomics.
</p>

### PanGenome Graph Evaluator [(PGGE)](https://github.com/pangenome/pgge)
----> Brief Description <----

---
# Learning Objectives

1. Creating scripts in specific work directory in the Nesi environment.
2. Downloading sequencing data (in fastq format). 
3. Creating Pangenome graphs using [PGGB](https://github.com/pangenome/pggb) and [PGGE](https://github.com/pangenome/pgge) Tools.
4. Identifying genetic variation and structural variance.

# Working in the Nesi Environment 
NeSI HPC environment is used for the analysis. Please make sure to have a NeSI account and you are able to login.

### Setting up your project directory

```bash
# Create a new directory in somewhere and change to that directory
mkdir pg_workshop
cd pg_workshop
# Keep a note of the absolute path of your directory
pwd
/nesi/nobackup/ga03793/pg_workshop
```

# Genome Data

### Genome Availability 
----> Mention the data source <----

### Other Data Availability
----> Mention other data source <----

---
# Procedure 
### 1. Downloading the Assembly Data file 4Sim.fa

```bash
#Use wget to download the file from the URL
wget https://github.com/ZoeYang2020/Pangenome-graph-for-bacterial-pathogens/raw/main/ESR_pangenome_pipeline_v2.0/4Sim.fa
```

### 2. Creating an index for the seuqence file and check

```bash
#Use samtools to create the index file
#In Nesi environment you will have to load the command first

module load SAMtools

samtools faidx 4Sim.fa 
ls -ltrh
total 8.8M
-rw-rw-r-- 1 ismnu81p ismnu81p 8.8M Apr 19 16:09 4Sim.fa
-rw-rw-r-- 1 ismnu81p ismnu81p  120 Apr 19 16:20 4Sim.fa.fai

cat 4Sim.fa.fai
NC_neisseria    2248966 14      60      61
Sim1_3k 2249014 2286472 60      61
Sim2_4k 2248965 4572979 60      61
Sim3_5k 2249050 6859436 60      61
```

As per the index this assembly consists of 4 samples described in the below table. 

| Name | Length |
|------|-------:|
|NC_neisseria | 2,248,966 |
|Sim1_3k | 2,249,014 |
|Sim2_4k | 2,248,965 |
|Sim3_5k | 2,249,050 |

### 3. Executing `pggb` tool using Singularity container
We can follow the procedure in https://github.com/pangenome/pggb#singularity to setup the Singularity image. This is already done and the image is in `/nesi/project/ga03793/software/pggb/` directory for version 0.5.3. 

Following script ([pggb_test.sh](https://github.com/nuzla/Pangenome-Graphs-Workshop/blob/main/Scripts/pggb_test.sh)) can be used to run `pggb` on the downloaded sequence. 

```sh
#!/usr/bin/bash
module load Singularity
#export container to a variable for convenience
WD=/nesi/nobackup/ga03793/pg_workshop #Working Directory
container=/nesi/project/ga03793/software/pggb/pggb_0.5.3.simg
data=${WD}/4Sim.fa

#Bind filesystem to container image 
export SINGULARITY_BIND="${WD}, /nesi/project/ga03793/"

singularity exec ${container} pggb -i $data -s 1000 -p 95 -n 8 -k 79 -t 2 -S -m -o output 
```

In `pggb` `-i` is for specifying the sequence file. `-s` specifies the segment length for mapping and `-p` specifies percent identity for mapping/alignment. `-n` is for number of haplotypes. `-k` for minimum matching length. `-t` says number of threads to be used for the execution. `-S` will generate the stats. `-m` will generate MultiQC report of graphs' statistics and visualizations. `-o` specifies the output directory name. 



<!---
### 2.	QC and Adapter Trimming
### 3.	Alignment to Reference 
### 4.	Reference SNPs
### 5.	Other Analysis 
--->


