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
<p align="justify">
The pangenome graph evaluation is a pipeline in the pangenome graphs, which measures the reconstruction accuracy of the graph. It helps find the best pangenome graph using input data and tasks.
</p>

#### Learning Objectives

1. Creating scripts in specific work directory in the Nesi environment.
2. Downloading and preparing sequencing data (in fasta format). 
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
The National Library of Medicine is the largest library focused on biomedicine worldwide, serving as the central hub for biomedical informatics and computational biology. It has many genome assembly data and [Genome assembly ASM19152v1](https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_000191525.1/) will be used for this workshop. 

---
# Procedure 
### 1. Downloading and preparng assembly data file 4Sim.fa

Please follow the proedure described in this [page](https://github.com/nuzla/Pangenome-Graphs-Workshop/blob/main/preparing_data_files.md)

### 2. Creating an index for the seuqence file and check

```bash
#Use samtools to create the index file
#In Nesi environment you will have to load the command first

$ module load SAMtools

$ samtools faidx ASM19152v1_pgsim.fa 

$ $ cat ASM19152v1_pgsim.fa.fai 
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
|NC_017518.1_SNP_5000                | 2,248,966 |   5,000|       0|   0 |   0 |
|NC_017518.1_INDEL_5000               | 2,249,048 |       0|   5,000|   0 |   0 |
|NC_017518.1_SNP_4000_INDEL_4000      | 2,153,883 |   4,000|   4,000|   0 |   0 |
|NC_017518.1_SNP_4000_INDEL_4000_INV_4| 2,242,147 |   4,000|   4,000|   4 |   0 |
|NC_017518.1_SNP_4000_INDEL_4000_CNV_4| 2,415,498 |   4,000|   4,000|   0 |   4 |

### 3. Executing `pggb` tool using Singularity container
We can follow the procedure in https://github.com/pangenome/pggb#singularity to setup the Singularity image. This is already done and the image is in `/nesi/project/ga03793/software/pggb/` directory for version 0.5.3. 

Following script ([pggb_test.sh](https://github.com/nuzla/Pangenome-Graphs-Workshop/blob/main/Scripts/pggb_test.sh)) can be used to run `pggb` on the downloaded sequence. 

```bash
#!/usr/bin/bash
module load Singularity
#export container to a variable for convenience
WD=/nesi/nobackup/ga03793/pg_workshop #Working Directory
container=/nesi/project/ga03793/software/pggb/pggb_0.5.3.simg
data=${WD}/ASM19152v1_pgsim.fa

#Bind filesystem to container image 
export SINGULARITY_BIND="${WD}, /nesi/project/ga03793/"

singularity exec ${container} pggb -i $data -s 1000 -p 95 -n 6 -k 79 -t 2 -S -m -o output -V 'NC_017518.1:#' 
```

In `pggb` `-i` is for specifying the sequence file. `-s` specifies the segment length for mapping and `-p` specifies percent identity for mapping/alignment. `-n` is for number of haplotypes (or number of samples). `-k` for minimum matching length. `-t` says number of threads to be used for the execution. `-S` will generate the stats. `-m` will generate MultiQC report of graphs' statistics and visualizations. `-o` specifies the output directory name. `-V 'NC_017518.1:#'` will create a vcf file and its stats considering NC_017518.1 as the reference sequence. 

---
# Executing `pggb` as a [SLURM](https://github.com/SchedMD/slurm) Job

Executing shell scripts in the Nesi environment might not be the best way to handle larger files which will require large memory, CPU power and time. We can modify the previusely explained script as below ([pggb_slurm_1K95.sh](https://github.com/nuzla/Pangenome-Graphs-Workshop/blob/main/Scripts/pggb_slurm_1K95.sh)) to run as SLURM job. Note the additional parameters specified by `#SBATCH` which will indicate maximum resource limitations. 

```bash
#!/usr/bin/bash

#SBATCH --account       ga03793
#SBATCH --job-name      NC_017518.1_1K95
#SBATCH --cpus-per-task 8 
#SBATCH --mem           4G
#SBATCH --time          1:00:00

module purge
module load Singularity

#export container to a variable for convenience
WD=/nesi/nobackup/ga03793/pg_workshop #Working Directory
container=/nesi/project/ga03793/software/pggb/pggb_0.5.3.simg
data=${WD}/ASM19152v1_pgsim.fa

#Bind filesystem to container image 
export SINGULARITY_BIND="${WD}, /nesi/project/ga03793/"

singularity exec ${container} pggb -i $data -s 1000 -p 95 -n 6 -k 79 -t 24 -S -m -o output_1K95 -V 'NC_017518.1:#'  
```

The job can be submitted using the `sbatch` command it will show a job id. In this case 35887085

```
$ sbatch pggb_slurm_1K95.sh
Submitted batch job 35887085
```

We can monitor the job status using `seff` and `squeue` specifying the job id. 

```
seff 35887085
Job ID: 35887085
Cluster: mahuika
User/Group: ismnu81p/ismnu81p
State: RUNNING
Nodes: 1
Cores per node: 8
CPU Utilized: 00:00:00
CPU Efficiency: 0.00% of 00:04:16 core-walltime
Job Wall-clock time: 00:00:32
Memory Utilized: 0.00 MB (estimated maximum)
Memory Efficiency: 0.00% of 4.00 GB (4.00 GB/node)
WARNING: Efficiency statistics may be misleading for RUNNING jobs.
```

```
$ squeue --job 35887085
JOBID         USER     ACCOUNT   NAME        CPUS MIN_MEM PARTITI START_TIME     TIME_LEFT STATE    NODELIST(REASON)    
35887085      ismnu81p ga03793   NC_017518.1_   8      4G large   2023-05-21T0       58:35 RUNNING  wbn063 
```

SLURM will also create a output log file and we can monitor it realtime using  `tail -f`. 

```bash
$ tail -f slurm-35887085.out 
[smoothxg::(1-3)::prep] writing graph output_1K95/ASM19152v1_pgsim.fa.2ab4142.c2fac19.seqwish.gfa.prep.0.gfa
[smoothxg::(1-3)::main] building xg index
[smoothxg::(1-3)::smoothable_blocks] computing blocks
[smoothxg::(1-3)::smoothable_blocks] computing blocks for 54747 handles: 100.00% @ 5.47e+04/s elapsed: 00:00:00:01 remain: 00:00:00:00
[smoothxg::(1-3)::break_and_split_blocks] cutting blocks that contain sequences longer than max-poa-length (1400) and depth >= 0
[smoothxg::(1-3)::break_and_split_blocks] splitting 3862 blocks at identity 0.950 (WFA-based clustering) and at estimated-identity 0.950 (mash-based clustering)
[smoothxg::(1-3)::break_and_split_blocks] cutting and splitting 3862 blocks: 100.00% @ 1.49e+04/s elapsed: 00:00:00:00 remain: 00:00:00:00
[smoothxg::(1-3)::break_and_split_blocks] cut 0 blocks of which 0 had repeats
[smoothxg::(1-3)::break_and_split_blocks] split 0 blocks
[smoothxg::(1-3)::smooth_and_lace] applying local SPOA to 3862 blocks: 78.40% @ 5.77e+01/s elapsed: 00:00:00:52 remain: 00:00:00:14
```

When the job is completed the `seff` command will show a summary report with below details. The job has used 785.61 MB memory and taken 7 minuted and 50 seconds to complete. 

```
$ seff 35887085
Job ID: 35887085
Cluster: mahuika
User/Group: ismnu81p/ismnu81p
State: COMPLETED (exit code 0)
Nodes: 1
Cores per node: 8
CPU Utilized: 00:41:23
CPU Efficiency: 66.04% of 01:02:40 core-walltime
Job Wall-clock time: 00:07:50
Memory Utilized: 785.61 MB
Memory Efficiency: 19.18% of 4.00 GB
```

Now we can try the same script by changing the `pggb` parameters `-s`, `-p` and `-k` and compare the results. Please refer [script folder](https://github.com/nuzla/Pangenome-Graphs-Workshop/tree/main/Scripts) for scripts. 

---
# MultiQC Report
The script generated [output](https://github.com/nuzla/Pangenome-Graphs-Workshop/tree/main/Output) directory consists of a compehensive and interactive [MutltiQC Report](https://multiqc.info/) which will decribe all. Open the file multiqc_report.html which is in the output folder from your browser.

_Note: To download the output folder from the Nesi environment you can first zip it using the command `zip -r output.zip output`_

## Graph Viszualization Details (`-s 1000`, `-p 95`)

### ODGI Compressed 1D visualization
![ODGI Compressed 1D visualization](https://github.com/nuzla/Pangenome-Graphs-Workshop/blob/main/Output/ASM19152v1_pgsim.fa.2ab4142.c2fac19.c47d9e7.smooth.final.og.viz_O_multiqc.png?raw=true)

### ODGI 1D visualization
![ODGI 1D visualization](https://github.com/nuzla/Pangenome-Graphs-Workshop/blob/main/Output/ASM19152v1_pgsim.fa.2ab4142.c2fac19.c47d9e7.smooth.final.og.viz_multiqc.png?raw=true)

### ODGI 1D visualization by path position
![ODGI 1D visualization by path position](https://github.com/nuzla/Pangenome-Graphs-Workshop/blob/main/Output/ASM19152v1_pgsim.fa.2ab4142.c2fac19.c47d9e7.smooth.final.og.viz_pos_multiqc.png?raw=true)

### ODGI 1D visualization by path orientation
![ODGI 1D visualization by path orientation](https://github.com/nuzla/Pangenome-Graphs-Workshop/blob/main/Output/ASM19152v1_pgsim.fa.2ab4142.c2fac19.c47d9e7.smooth.final.og.viz_inv_multiqc.png?raw=true)

### ODGI 1D visualization by node depth
![ODGI 1D visualization by node depth](https://github.com/nuzla/Pangenome-Graphs-Workshop/blob/main/Output/ASM19152v1_pgsim.fa.2ab4142.c2fac19.c47d9e7.smooth.final.og.viz_depth_multiqc.png?raw=true)

### ODGI 1D visualization by uncalled bases
![ODGI 1D visualization by uncalled bases](https://github.com/nuzla/Pangenome-Graphs-Workshop/blob/main/Output/ASM19152v1_pgsim.fa.2ab4142.c2fac19.c47d9e7.smooth.final.og.viz_uncalled_multiqc.png?raw=true)

### ODGI 2D drawing
![ODGI 2D drawing](https://github.com/nuzla/Pangenome-Graphs-Workshop/blob/main/Output/ASM19152v1_pgsim.fa.2ab4142.c2fac19.c47d9e7.smooth.final.og.lay.draw.png?raw=true)

# Variant Call Analysis using the VCF file (`-s 1000`)

### Finding stats of the VCF file
`bcftools stats <file.vcf>` command will display the all the stats related to the VCF file. 

### Variant Call Comaprison

#### 1. Creating VCF file using linear method
The procedure described in [this page](https://github.com/nuzla/Pangenome-Graphs-Workshop/blob/main/linear_reference_vc.md) can be used to find linear reference based stats. 

#### 2. Creating VCF file using graph method
_**This section work in progress**_

As explained in a previouse section we sepcified the option `-V 'NC_017518.1:#'`. That will execute the command,
```
vg deconstruct -P NC_017518.1 -H # -e -a -t 2 output/ASM19152v1_pgsim.fa.2ab4142.c2fac19.c47d9e7.smooth.final.gfa
```
You can see it in the [Log File](https://github.com/nuzla/Pangenome-Graphs-Workshop/blob/main/Output/ASM19152v1_pgsim.fa.2ab4142.c2fac19.c47d9e7.smooth.05-21-2023_04_21_52.log) line 279.

### Sample wise stats

_**This section work in progress**_

#### 1. NC_017518.1_SNP_5000 (2,248,966)

|Type	        |    SNP	    |  INDEL	 |    INV	  |    CNV	 |   Total	  |    TP	   |    TN	   |    FP	   |     FN	   |  Sensitivity	| Specificity |
|---------:   |---------:  |--------:|---------:|--------:|---------: |---------:|---------:| --------:|--------:  |---------:    |--------:    |
|Ground Truth |    5,000	  |  0    	 |    0  	  |    0  	 |   5,000	  |      	   |      	   |      	   |      	    |   	          |             |
|Linear Method|    5,000	  |  0    	 |    0  	  |    0  	 |   5,000	  |      	   |      	   |      	   |      	    |   	          |             |
|pggb (95%)   |         	  |       	 |       	  |       	 |           |      	   |      	   |      	   |      	    |   	          |             |
|pggb (98%)   |         	  |       	 |       	  |       	 |           |      	   |      	   |      	   |      	    |   	          |             |

#### 2. NC_017518.1_INDEL_5000	(2,249,048)


|Type	        |    SNP	    |  INDEL	 |    INV	  |    CNV	 |   Total	  |    TP	   |    TN	   |    FP	   |     FN	   |  Sensitivity	| Specificity |
|---------:   |---------:  |--------:|---------:|--------:|---------: |---------:|---------:| --------:|--------:  |---------:    |--------:    |
|Ground Truth |       0 	  |  5,000	 |    0  	  |    0  	 |   5,000	  |      	   |      	   |      	   |      	    |   	          |             |
|Linear Method|    329  	  |  3,977 	|    0  	  |    0  	 |   4,306	  |      	   |      	   |      	   |      	    |   	          |             |
|pggb (95%)   |         	  |       	 |       	  |       	 |           |      	   |      	   |      	   |      	    |   	          |             |
|pggb (98%)   |         	  |       	 |       	  |       	 |           |      	   |      	   |      	   |      	    |   	          |             |

#### 3. NC_017518.1_SNP_4000_INDEL_4000 (2,153,883)

|Type	        |    SNP	    |  INDEL	 |    INV	  |    CNV	 |   Total	  |    TP	   |    TN	   |    FP	   |     FN	   |  Sensitivity	| Specificity |
|---------:   |---------:  |--------:|---------:|--------:|---------: |---------:|---------:| --------:|--------:  |---------:    |--------:    |
|Ground Truth |    4,000	  |  4,000  |    0  	  |    0  	 |   8,000	  |      	   |      	   |      	   |      	    |   	          |             |
|Linear Method|    4,022	  |  3,282  |    0  	  |    0  	 |   7,304	  |      	   |      	   |      	   |      	    |   	          |             |
|pggb (95%)   |         	  |       	 |       	  |       	 |           |      	   |      	   |      	   |      	    |   	          |             |
|pggb (98%)   |         	  |       	 |       	  |       	 |           |      	   |      	   |      	   |      	    |   	          |             |

#### 4. NC_017518.1_SNP_4000_INDEL_4000_INV_4	(2,242,147)

|Type	        |    SNP	    |  INDEL	 |    INV	  |    CNV	 |   Total	  |    TP	   |    TN	   |    FP	   |     FN	   |  Sensitivity	| Specificity |
|---------:   |---------:  |--------:|---------:|--------:|---------: |---------:|---------:| --------:|--------:  |---------:    |--------:    |
|Ground Truth |    4,000	  |  4,000  |    4  	  |    0  	 |   8,000	  |      	   |      	   |      	   |      	    |   	          |             |
|Linear Method|    4,021	  |  3,282  |    0  	  |    0  	 |   7,303	  |      	   |      	   |      	   |      	    |   	          |             |
|pggb (95%)   |         	  |       	 |       	  |       	 |           |      	   |      	   |      	   |      	    |   	          |             |
|pggb (98%)   |         	  |       	 |       	  |       	 |           |      	   |      	   |      	   |      	    |   	          |             |

#### 5. NC_017518.1_SNP_4000_INDEL_4000_CNV_4	(2,415,498)

|Type	        |    SNP	    |  INDEL	 |    INV	  |    CNV	 |   Total	  |    TP	   |    TN	   |    FP	   |     FN	   |  Sensitivity	| Specificity |
|---------:   |---------:  |--------:|---------:|--------:|---------: |---------:|---------:| --------:|--------:  |---------:    |--------:    |
|Ground Truth |    4,000	  |  4,000  |    0  	  |    4  	 |   8,000	  |      	   |      	   |      	   |      	    |   	          |             |
|Linear Method|    4,111   |  3,247  |    0  	  |    0  	 |   7,358	  |      	   |      	   |      	   |      	    |   	          |             |
|pggb (95%)   |         	  |       	 |       	  |       	 |           |      	   |      	   |      	   |      	    |   	          |             |
|pggb (98%)   |         	  |       	 |       	  |       	 |           |      	   |      	   |      	   |      	    |   	          |             |


<!---
#### 1. Linear Reference VCF vs PGGB Graph VCF.

The procedure described in [this page](https://github.com/nuzla/Pangenome-Graphs-Workshop/blob/main/linear_reference_vc.md) can be used to find linear reference based stats and `bcftools isec` can be used to compare VCF files. 

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

\end{aligned}
```
Now repeat the procedure for Identity 90% and 95%. 

|Identity| Linear Ref SNP   | PGGB Graph SNP | True Positive (TP) | False Positive (FP) | True Negative (TN) | False Negative (FN) | Sensitivity | Specifisity|
|-------:|-----------------:|---------------:|-------------------:|--------------------:|-------------------:|--------------------:|------------:|-----------:|
|90%     |         10,041   |           8,606|               8,546|                   60|                    |                1,495|       85.11%|            |
|95%     |         10,041   |           8,678|               8,617|                   61|                    |                1,424|       85.81%|            |
|98%     |         10,041   |           8,617|               8,553|                   64|                    |                1,488|       85.18%|            |
-->
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
