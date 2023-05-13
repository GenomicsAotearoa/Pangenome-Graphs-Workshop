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

singularity exec ${container} pggb -i $data -s 1000 -p 95 -n 4 -k 79 -t 2 -S -m -o output -V 'NC_neisseria:#'
```

In `pggb` `-i` is for specifying the sequence file. `-s` specifies the segment length for mapping and `-p` specifies percent identity for mapping/alignment. `-n` is for number of haplotypes (or number of samples). `-k` for minimum matching length. `-t` says number of threads to be used for the execution. `-S` will generate the stats. `-m` will generate MultiQC report of graphs' statistics and visualizations. `-o` specifies the output directory name. `-V 'NC_neisseria:#'` will create a vcf file and it stats considering NC_neisseria as the reference sequence. 

---
# Executing `pggb` as a [SLURM](https://github.com/SchedMD/slurm) Job

Executing shell scripts in the Nesi environment might not be the best way to handle larger files which will require large memory, CPU power and time. We can modify the previusely explained script as below ([pggb_slurm_1K95.sh](https://github.com/nuzla/Pangenome-Graphs-Workshop/blob/main/Scripts/pggb_slurm_1K95.sh)) to run as SLURM job. Note the additional parameters specified by `#SBATCH` which will indicate maximum resource limitations. 

```bash
#!/usr/bin/bash

#SBATCH --account       ga03793
#SBATCH --job-name      4Sim_1K95
#SBATCH --cpus-per-task 4 
#SBATCH --mem           4G
#SBATCH --time          1:00:00

module purge
module load Singularity

#export container to a variable for convenience
WD=/nesi/nobackup/ga03793/pg_workshop #Working Directory
container=/nesi/project/ga03793/software/pggb/pggb_0.5.3.simg
data=${WD}/4Sim.fa

#Bind filesystem to container image 
export SINGULARITY_BIND="${WD}, /nesi/project/ga03793/"

singularity exec ${container} pggb -i $data -s 1000 -p 95 -n 4 -k 79 -t 2 -S -m -o 4Sim_1K95 -V 'NC_neisseria:#' 
```

The job can be submitted using the `sbatch` command it will show a job id. In this case 34588496

```bash
sbatch pggb_slurm_1K95.sh 
Submitted batch job 34588496
```

We can monitor the job status using `seff` and `squeue` specifying the job id. 

```bash
seff 34588496
Job ID: 34588496
Cluster: mahuika
User/Group: ismnu81p/ismnu81p
State: RUNNING
Nodes: 1
Cores per node: 4
CPU Utilized: 00:00:00
CPU Efficiency: 0.00% of 00:04:56 core-walltime
Job Wall-clock time: 00:01:14
Memory Utilized: 0.00 MB (estimated maximum)
Memory Efficiency: 0.00% of 4.00 GB (4.00 GB/node)
WARNING: Efficiency statistics may be misleading for RUNNING jobs.
```

```bash
squeue --job 34588496
JOBID         USER     ACCOUNT   NAME        CPUS MIN_MEM PARTITI START_TIME     TIME_LEFT STATE    NODELIST(REASON)    
34588496      ismnu81p ga03793   4Sim_1K95      4      4G large   2023-04-20T0       58:10 RUNNING  wbn182
```

SLURM will also create a output log file and we can monitor it realtime using  `tail -f`. 

```bash
tail -f slurm-34588496.out
[smoothxg::(1-3)::prep] writing graph 4Sim_1K95/4Sim.fa.3541aba.c2fac19.seqwish.gfa.prep.0.gfa
[smoothxg::(1-3)::main] building xg index
[smoothxg::(1-3)::smoothable_blocks] computing blocks
[smoothxg::(1-3)::smoothable_blocks] computing blocks for 36095 handles: 100.00% @ 7.21e+04/s elapsed: 00:00:00:00 remain: 00:00:00:00
[smoothxg::(1-3)::break_and_split_blocks] cutting blocks that contain sequences longer than max-poa-length (1400) and depth >= 0
[smoothxg::(1-3)::break_and_split_blocks] splitting 3459 blocks at identity 0.950 (WFA-based clustering) and at estimated-identity 0.950 (mash-based clustering)
[smoothxg::(1-3)::break_and_split_blocks] cutting and splitting 3459 blocks: 100.00% @ 1.37e+04/s elapsed: 00:00:00:00 remain: 00:00:00:00
[smoothxg::(1-3)::break_and_split_blocks] cut 0 blocks of which 0 had repeats
[smoothxg::(1-3)::break_and_split_blocks] split 0 blocks
[smoothxg::(1-3)::smooth_and_lace] applying local SPOA to 3459 blocks: 21.97% @ 2.91e+01/s elapsed: 00:00:00:26 remain: 00:00:01:32
```

When the job is completed the `seff` command will show a summary report with below details. The job has used 147.44 MB memory and taken 8 minuted and 25 seconds to complete. 

```bash
seff 34588496
Job ID: 34588496
Cluster: mahuika
User/Group: ismnu81p/ismnu81p
State: COMPLETED (exit code 0)
Nodes: 1
Cores per node: 4
CPU Utilized: 00:15:20
CPU Efficiency: 45.54% of 00:33:40 core-walltime
Job Wall-clock time: 00:08:25
Memory Utilized: 147.44 MB
Memory Efficiency: 3.60% of 4.00 GB
```

Now we can try the same script by changing the `pggb` parameters `-s`, `-p` and `-k` and compare the results. 

---
# MultiQC Report
The script generated [output](https://github.com/nuzla/Pangenome-Graphs-Workshop/tree/main/Output) directory consists of a compehensive and interactive [MutltiQC Report](https://multiqc.info/) which will decribe all. Open the file multiqc_report.html which is in the output folder from your browser.

_Note: To download the output folder from the Nesi environment you can first zip it using the command `zip -r output.zip output`_

## Graph Viszualization Details (`-s 1000`, `-p 95`)

### ODGI Compressed 1D visualization
![ODGI Compressed 1D visualization](https://github.com/nuzla/Pangenome-Graphs-Workshop/blob/main/Output/4Sim.fa.3541aba.c2fac19.9d98660.smooth.final.og.viz_O_multiqc.png?raw=true)

### ODGI 1D visualization
![ODGI 1D visualization](https://github.com/nuzla/Pangenome-Graphs-Workshop/blob/main/Output/4Sim.fa.3541aba.c2fac19.9d98660.smooth.final.og.viz_multiqc.png?raw=true)

### ODGI 1D visualization by path position
![ODGI 1D visualization by path position](https://github.com/nuzla/Pangenome-Graphs-Workshop/blob/main/Output/4Sim.fa.3541aba.c2fac19.9d98660.smooth.final.og.viz_pos_multiqc.png?raw=true)

### ODGI 1D visualization by path orientation
![ODGI 1D visualization by path orientation](https://github.com/nuzla/Pangenome-Graphs-Workshop/blob/main/Output/4Sim.fa.3541aba.c2fac19.9d98660.smooth.final.og.viz_inv_multiqc.png?raw=true)

### ODGI 1D visualization by node depth
![ODGI 1D visualization by node depth](https://github.com/nuzla/Pangenome-Graphs-Workshop/blob/main/Output/4Sim.fa.3541aba.c2fac19.9d98660.smooth.final.og.viz_depth_multiqc.png?raw=true)

### ODGI 1D visualization by uncalled bases
![ODGI 1D visualization by uncalled bases](https://github.com/nuzla/Pangenome-Graphs-Workshop/blob/main/Output/4Sim.fa.3541aba.c2fac19.9d98660.smooth.final.og.viz_uncalled_multiqc.png?raw=true)

### ODGI 2D drawing
![ODGI 2D drawing](https://github.com/nuzla/Pangenome-Graphs-Workshop/blob/main/Output/4Sim.fa.3541aba.c2fac19.9d98660.smooth.final.og.lay.draw.png?raw=true)

# SNP Analysis using the VCF file (`-s 1000`)

### How `pggb` creates VCF file?
As explained in a previuse section we sepcified the option `-V 'NC_neisseria:#'`. That will execute the command,
```
vg deconstruct -P NC_neisseria -H # -e -a -t 2 output/4Sim.fa.3541aba.c2fac19.9d98660.smooth.final.gfa
```
You can see it in the [Log File](https://github.com/nuzla/Pangenome-Graphs-Workshop/blob/main/Output/4Sim.fa.3541aba.c2fac19.9d98660.smooth.04-20-2023_02_56_38.log) line 296. 

### How to create the Ground Truth VCF?
This procedure is exaplned in [this page](https://github.com/nuzla/Pangenome-Graphs-Workshop/blob/main/ground_truth_data.md).

### Finding stats of the VCF file
`bcftools stats <file.vcf>` command will display the all the stats related to the VCF file. 

### SNP Comaprison

#### 1. Ground truth VCF vs PGGB VCF.

In order to find the True Positive (TP), True Negative (TN), False Positive (FP) and False Negative (FN) we need to compare the 2 VCF files. `bcftools isec` can be used for that. 
e.g. :
```
bcftools isec -p output_dir 4Sim_ref.vcf.gz 4Sim_1K98.vcf.gz 
cd output_dir/
cat README.txt 
This file was produced by vcfisec.
The command line was:	bcftools isec  -p output_dir 4Sim_ref.vcf.gz 4Sim_1K98.vcf.gz

Using the following file names:
output_dir/0000.vcf	for records private to	4Sim_ref.vcf.gz
output_dir/0001.vcf	for records private to	4Sim_1K98.vcf.gz
output_dir/0002.vcf	for records from 4Sim_ref.vcf.gz shared by both	4Sim_ref.vcf.gz 4Sim_1K98.vcf.gz
output_dir/0003.vcf	for records from 4Sim_1K98.vcf.gz shared by both	4Sim_ref.vcf.gz 4Sim_1K98.vcf.gz
```
As usual we can apply `bcftools stats <file.vcf>` to get the stats form each file. We can illustrate it with the below Venn diagram.

<img src="https://github.com/nuzla/Pangenome-Graphs-Workshop/assets/8539123/e9cfe33b-e473-4a25-a190-7af2e603b654" width="600">

For this case False Negative (Undetected in Graph) count is 1,488 and False Postive (Wrongly Detected in Graph) count is 64. 

```math
\begin{aligned}
Sensitivity  & = \frac{TP}{TP+FN} \\
              &  = \frac{8553}{8553+1488} \\
              & = 85.18\% \\ \\
Specificity & = \frac{TN}{TN+FP} \\
            & = \frac{496}{496+64} \\
            & = 88.57\%
\end{aligned}
```

TN considered as (number of records - number of SNPs) of record from reference and shared by both. 

Now repeat the procedure for Identity 90% and 95%. 

|Identity| Ground Truth SNP | PGGB Graph SNP | True Positive (TP) | False Positive (FP) | True Negative (TN) | False Negative (FN) | Sensitivity | Specifisity|
|-------:|-----------------:|---------------:|-------------------:|--------------------:|-------------------:|--------------------:|------------:|-----------:|
|90%     |         10,041   |           8,606|               8,546|                   60|                 512|                1,495|       85.11%|      89.51%|
|95%     |         10,041   |           8,678|               8,617|                   61|                 498|                1,424|       85.81%|      89.08%|
|98%     |         10,041   |           8,617|               8,553|                   64|                 496|                1,488|       85.18%|      88.57%|

# PanGenome Graph Evaluator
---
## Executing pgge tool using Singularity container
First we need to download a .R file to the working directory. 
```
wget https://raw.githubusercontent.com/pangenome/pgge/master/scripts/beehave.R
```
Now execute the script [pgge_test.sh](https://github.com/nuzla/Pangenome-Graphs-Workshop/blob/main/Scripts/pgge_test.sh)

```bash
#!/usr/bin/bash
module load Singularity
#export container to a variable for convenience
WD=/nesi/nobackup/ga03793/pg_workshop #Working Directory
container=/nesi/project/ga03793/software/pgge/pgge_032023.simg
data=${WD}/4Sim.fa

#Bind filesystem to container image 
export SINGULARITY_BIND="${WD}, /nesi/project/ga03793/"

singularity exec ${container} pgge -g ${WD}/output/*.gfa -f $data -o pgge -r ${WD}/beehave.R -b pgge/pgge_4Sim_peanut_bed -l 100000 -s 5000 -t 16
```
---> Please explain the options <----

You can find the result out pgge folder [here](https://github.com/nuzla/Pangenome-Graphs-Workshop/tree/main/pgge)

You can also run this task as a SLURM job using the script [pgge_slurm.sh](https://github.com/nuzla/Pangenome-Graphs-Workshop/blob/main/Scripts/pgge_slurm.sh)

The output graph from `pgge` is shown below. 

![pgge-l100000-s5000.pgge.tsv.png](https://github.com/nuzla/Pangenome-Graphs-Workshop/blob/main/pgge/pgge-l100000-s5000.pgge.tsv.png?raw=true)
