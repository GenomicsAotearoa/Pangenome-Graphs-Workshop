# Pangenome Graphs Overview
This repository contains material to construct pangenome graphs for a small bacteria dataset detailing every step in the Nesi environment. 

A pangenome graph is a graphical representation of the collective genomic information of a set of related organisms. Unlike a traditional linear genome assembly, which represents a single consensus genome sequence, a pangenome graph captures the genetic variation and structural diversity within a population. Pangenome graphs are constructed by integrating multiple genome sequences into a single graph structure that represents the entire set of genetic elements present in the population. This graph structure allows for the identification and visualization of genomic variations, such as single nucleotide polymorphisms (SNPs), insertions, deletions, and structural variations, as well as the relationships between different genomic regions. Pangenome graphs have become an important tool in genomics research, especially in the study of bacterial and viral populations, where genetic diversity is often high.

Neisseria is a genus of Gram-negative bacteria that are typically found in the mucous membranes and epithelial tissues of humans and animals. There are various species of Neisseria bacteria, some of which are harmless commensals, while others can cause diseases. Two well-known pathogenic species are Neisseria meningitidis, which can cause meningitis and septicemia, and Neisseria gonorrhoeae, which is the causative agent of the sexually transmitted infection gonorrhea. Other species of Neisseria are generally considered harmless and reside as commensals in the oral and/or nasopharynx of humans.

To investigate and analyse Neisseria Bacteria the PanGenome Graph Builder has been used. PanGenome Graph Builder (PGGB) is a computational tool used in genomics research to construct pan-genome graphs from large sets of genomic data. A pan-genome is the collection of all the genes and non-coding sequences present in a species or a group of related organisms. PGGB constructs a pan-genome graph, which is a data structure that represents the entire set of genes and genetic variations in a population. This graph can be used to study genetic diversity, gene function, and evolution. PGGB is designed to be scalable and efficient, making it suitable for large-scale genomic analyses. It is an open-source tool that can be used freely by researchers in the field of genomics.

This study includes an anaysis of Neisseria Bacteria genome sequence data with 4Sim data samples to construct pangenome graphs to identify genetic variation and structural variance.

### Learning Objectives
1. Creating scripts in specific work directory in the Nesi environment.
2. Downloading sequencing data (in fastq format). 
3. Creating Pangenome graphs using [PGGB](https://github.com/pangenome/pggb){:target="_blank"} and [PGGE](https://github.com/pangenome/pgge){:target="_blank"} Tools.
4. Identifying genetic variation and structural variance.


# Working in the Nesi Environment 
NeSI HPC environment is used for the analysis. Please make sure to have a NeSI account and you are able to login.

### Setting up your project directories

```
# Create a directory to do the analysis, and change to that directory
mkdir Workshop
cd Workshop
```

# Genome Data
### Downloading the required sequence file. 
```
#Use wget to download the file from the URL
wget https://github.com/ZoeYang2020/Pangenome-graph-for-bacterial-pathogens/raw/main/ESR_pangenome_pipeline_v2.0/4Sim.fa
```

### Creating an index for the seuqence file and check
```
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

### Genome Availability 

### Other Data Availability
# Methods 
### 1.	Assembly Procedure 
### 2.	QC and Adapter Trimming
### 3.	Alignment to Reference 
### 4.	Reference SNPs
### 5.	Other Analysis 


# Scripts 

# Dependencies 
