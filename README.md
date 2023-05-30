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

### How does the pggb graph build work?
The pggb algorithm facilitates the construction of these graphs by progressively integrating genetic variants into a reference genome.
The key features and purpose of the PGGb include:
1.	Progressive Approach: The PGGb algorithm follows a step-by-step approach, iteratively adding genetic variants to the graph one at a time.
2.	Variant Integration: The tool efficiently integrates genetic variants, including single nucleotide variations (SNVs), insertions, deletions, and larger structural variations, into the graph representation.
3.	Optimal Placement: For each variant, the PGGb algorithm determines the most suitable position to incorporate it into the graph. This involves aligning the variant sequence with the existing graph structure while minimizing conflicts with the nodes and edges.
4.	Graph Expansion: Once the optimal placement is determined, the PGGb algorithm expands the graph by adding new nodes and edges to represent the variant sequence. The overall graph structure is modified to connect the variant nodes with the adjacent reference nodes.
5.	Large-Scale Graph Construction: The PGGb algorithm is designed to handle large-scale genomes and can efficiently construct genome graphs containing extensive genetic variations.

### All-to-all alignment
Generally, refers to the process of aligning all sequences in a given set against each other, rather than aligning them to a single reference sequence.
We begin with an alignment, with wfmash. This compares all sequences to each other and finds the best N mappings for each. It produces base-level alignments.

### Inducing the graph
Refers to the process of constructing the genome graph by progressively integrating genetic variants into a reference genome.
These base-level alignments are converted into a graph with seqwish. A filter is applied to remove short matches, which anchors the graph on confident longer exact matches.

### Normalizing the graph
refers to a process that aims to optimize the structure and representation of the genome graph by resolving redundant or overlapping elements. This step is typically performed after the initial construction of the graph.
The normalization process in PGGb involves several steps, which may vary depending on the specific implementation or version of the tool. Here are some common steps involved in normalizing the graph:
1.	Removal of Redundant Nodes: During the construction of the genome graph, it is possible that some nodes become redundant due to overlapping or repetitive sequences. Normalization involves identifying and removing these redundant nodes, streamlining the graph structure.
2.	Edge Optimization: In the graph, edges represent connections between nodes. During normalization, the edges are optimized to minimize redundancy and improve the efficiency of the graph. This can include merging or repositioning edges to create a more streamlined and accurate representation of the genome.
3.	Compact Representation: Normalization aims to reduce the overall size of the graph by compacting the representation. This can involve compressing repetitive regions or simplifying complex structures while preserving the essential information and variant representation.
4.	Graph Refinement: The normalization process also involves refining the graph structure by resolving inconsistencies, correcting errors, and improving the overall quality of the graph representation. This may include resolving conflicts between nodes and edges, addressing mismatches, and ensuring the graph accurately reflects the underlying genetic variations

To normalize the graph and harmonize the allele representation, we use smoothxg to apply a local MSA across all parts of the graph.

### Downstream
These graphs offer a wide range of capabilities. Initially, we can generate several diagnostic visualizations derived from the graphs, providing a user-friendly way to comprehend the alignments at a broader level. Additionally, using the PGGb tool, we can generate variant calls by leveraging vg deconstruct. The graphs produced by PGGb serve as reference systems for aligning short reads through vg giraffe or long reads through GraphAligner. Furthermore, odgi allows us to utilize the graphs as reference systems to elucidate homology relationships across entire genomes.

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
### 1. Downloading and preparing assembly data file 4Sim.fa

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

You can run run pggb without parameters to get information on the meaning of each parameter. Noe take a look at the files in the "output" folder.
We get a graph in GFA (*.gfa) and odgi (*.og) formats. These can be used downstream in many methods, including those in vg, like vg giraffe. You can visualize the GFA format graph with BandageNG, and use odgi directly on the \*.gfa or \*.og output. 

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

#### Understanding odgi visualizations
We obtain a series of diagnostic images that represent the pangenome alignment. These are created with odgi viz (1D matrix) and odgi layout with odgi draw (2D graph drawings). First, the 2D layout gives us a view of the total alignment. For small graphs, we can look at the version that shows where specific paths go (\*.draw_multiqc.png): For larger ones, the \*.draw.png result is usually more legible.

#### The effect of haplotype count `-n`
What happens if we set a lower `-n`? This parameter determines how many mappings we have. Each sequence is aligned against its n-1 best matches. Setting -n 6 causes clustering of sequences into groups that are more similar.

#### The effect of the minimum match filter `-k`
Another key parameter is -k, which affects the behavior of seqwish. This filter removes exact matches from alignments that are shorter than -k. Short matches occur in regions of high diversity. In practice, these short matches contribute little to the overall structure of the graph, and we remove them to further simplify the base graph structure.

### Decreasing mapping segment length `-s` increases sensitivity
By setting a lower mapping segment length, which affects the behavior of the very first step in the pipeline, wfmash’s mapping step (itself based on a heavily modified version of MashMap). This defaults to -s 5k. We can use -s 1k to guarantee we pick up on smaller homology segments, leading to a more complete alignment.

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












<!---
# Variant Call Analysis using the VCF file (`-s 1000`)

### Finding stats of the VCF file
`bcftools stats <file.vcf>` command will display the all the stats related to the VCF file. 

### Variant Call Comaprison

#### 1. Creating VCF file using linear method
The procedure described in [this page](https://github.com/nuzla/Pangenome-Graphs-Workshop/blob/main/linear_reference_vc.md) can be used to find linear reference based stats. 


#### 2. Creating VCF file using graph method
As explained in a previouse section we sepcified the option `-V 'NC_017518.1:#'`. That will execute the command,
```
vg deconstruct -P NC_017518.1 -H # -e -a -t 2 output/ASM19152v1_pgsim.fa.2ab4142.c2fac19.c47d9e7.smooth.final.gfa
```
You can see it in the [Log File](https://github.com/nuzla/Pangenome-Graphs-Workshop/blob/main/Output/ASM19152v1_pgsim.fa.2ab4142.c2fac19.c47d9e7.smooth.05-21-2023_04_21_52.log) line 279.

The generated VCF file is availble in the output folder. 

```
$ ls -1sh output/*.vcf
5.3M output/ASM19152v1_pgsim.fa.2ab4142.c2fac19.c47d9e7.smooth.final.NC_017518.1.vcf

$ less output/ASM19152v1_pgsim.fa.2ab4142.c2fac19.c47d9e7.smooth.final.NC_017518.1.vcf
##fileformat=VCFv4.2
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##INFO=<ID=AC,Number=A,Type=Integer,Description="Total number of alternate alleles in called genotypes">
##INFO=<ID=AF,Number=A,Type=Float,Description="Estimated allele frequency in the range (0,1]">
##INFO=<ID=NS,Number=1,Type=Integer,Description="Number of samples with data">
##INFO=<ID=AN,Number=1,Type=Integer,Description="Total number of alleles in called genotypes">
##INFO=<ID=LV,Number=1,Type=Integer,Description="Level in the snarl tree (0=top level)">
##INFO=<ID=PS,Number=1,Type=String,Description="ID of variant corresponding to parent snarl">
##INFO=<ID=AT,Number=R,Type=String,Description="Allele Traversal as path in graph">
##contig=<ID=NC_017518.1_INDEL_5000,length=2249048>
##contig=<ID=NC_017518.1_SNP_4000_INDEL_4000_INV_4,length=2242147>
##contig=<ID=NC_017518.1,length=2248966>
##contig=<ID=NC_017518.1_SNP_4000_INDEL_4000_CNV_4,length=2415498>
##contig=<ID=NC_017518.1_SNP_5000,length=2248966>
##contig=<ID=NC_017518.1_SNP_4000_INDEL_4000,length=2242147>
#CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT
NC_017518.1     195     >1>4    T       G       60      .       AC=0;AF=0;AN=0;AT=>1>3>4,>1>2>4;NS=0;LV=0       GT
NC_017518.1     295     >4>6    G       GATCTACCCTGCTA  60      .       AC=0;AF=0;AN=0;AT=>4>6,>4>5>6;NS=0;LV=0 GT
NC_017518.1     315     >6>9    C       A       60      .       AC=0;AF=0;AN=0;AT=>6>8>9,>6>7>9;NS=0;LV=0       GT
NC_017518.1     737     >9>11   A       AG      60      .       AC=0;AF=0;AN=0;AT=>9>11,>9>10>11;NS=0;LV=0      GT
NC_017518.1     814     >11>14  A       T       60      .       AC=0;AF=0;AN=0;AT=>11>13>14,>11>12>14;NS=0;LV=0 GT
NC_017518.1     872     >14>16  AG      A       60      .       AC=0;AF=0;AN=0;AT=>14>15>16,>14>16;NS=0;LV=0    GT
NC_017518.1     965     >16>19  A       C       60      .       AC=0;AF=0;AN=0;AT=>16>17>19,>16>18>19;NS=0;LV=0 GT
NC_017518.1     981     >19>21  TGG     T       60      .       AC=0;AF=0;AN=0;AT=>19>20>21,>19>21;NS=0;LV=0    GT
NC_017518.1     1485    >21>24  C       T       60      .       AC=0;AF=0;AN=0;AT=>21>23>24,>21>22>24;NS=0;LV=0 GT
```

But for our analysis we need path wise information belong to each sample. So we will rerun the `vg deconstruct` with some different parameters. 

```
$ vg deconstruct -p NC_017518.1 -a -e -K output/ASM19152v1_pgsim.fa.2ab4142.c2fac19.c47d9e7.smooth.final.gfa > ASM19152v1_pgsim.fa.smooth.final.NC_017518.1.vcf

$ less ASM19152v1_pgsim.fa.smooth.final.NC_017518.1.vcf
##fileformat=VCFv4.2
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##INFO=<ID=AC,Number=A,Type=Integer,Description="Total number of alternate alleles in called genotypes">
##INFO=<ID=AF,Number=A,Type=Float,Description="Estimated allele frequency in the range (0,1]">
##INFO=<ID=NS,Number=1,Type=Integer,Description="Number of samples with data">
##INFO=<ID=AN,Number=1,Type=Integer,Description="Total number of alleles in called genotypes">
##INFO=<ID=LV,Number=1,Type=Integer,Description="Level in the snarl tree (0=top level)">
##INFO=<ID=PS,Number=1,Type=String,Description="ID of variant corresponding to parent snarl">
##INFO=<ID=AT,Number=R,Type=String,Description="Allele Traversal as path in graph">
##contig=<ID=NC_017518.1,length=2248966>
#CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT  NC_017518.1_INDEL_5000  NC_017518.1_SNP_4000_INDEL_4000 NC_017518.1_SNP_4000_INDEL_4000_CNV_4   NC_017518.1_SNP_4000_INDEL_4000_INV_4   NC_017518.1_SNP_5000
NC_017518.1     195     >1>4    T       G       60      .       AC=3;AF=0.6;AN=5;AT=>1>3>4,>1>2>4;NS=5;LV=0     GT      0       1       1       1       0
NC_017518.1     295     >4>6    G       GATCTACCCTGCTA  60      .       AC=1;AF=0.2;AN=5;AT=>4>6,>4>5>6;NS=5;LV=0       GT      1       0       0       0       0
NC_017518.1     315     >6>9    C       A       60      .       AC=1;AF=0.2;AN=5;AT=>6>8>9,>6>7>9;NS=5;LV=0     GT      0       0       0       0       1
NC_017518.1     737     >9>11   A       AG      60      .       AC=1;AF=0.2;AN=5;AT=>9>11,>9>10>11;NS=5;LV=0    GT      1       0       0       0       0
NC_017518.1     814     >11>14  A       T       60      .       AC=3;AF=0.6;AN=5;AT=>11>13>14,>11>12>14;NS=5;LV=0       GT      0       1       1       1       0
NC_017518.1     872     >14>16  AG      A       60      .       AC=3;AF=0.6;AN=5;AT=>14>15>16,>14>16;NS=5;LV=0  GT      0       1       1       1       0
NC_017518.1     965     >16>19  A       C       60      .       AC=1;AF=0.2;AN=5;AT=>16>17>19,>16>18>19;NS=5;LV=0       GT      0       0       0       0       1
```
For viewing stats belongs to each path or sample we have to use some sample filter and non-variant row filter using the option `--min-ac=1`. To get the stats about the sample NC_017518.1_SNP_5000 only we have use the below syntax.

```
$ bcftools view -s NC_017518.1_SNP_5000 --min-ac=1 ASM19152v1_pgsim.fa.smooth.final.NC_017518.1.vcf | bcftools stats - | less 
# This file was produced by bcftools stats (1.9+htslib-1.9) and can be plotted using plot-vcfstats.
# The command line was: bcftools stats  -
#
# Definition of sets:
# ID    [2]id   [3]tab-separated file names
ID      0       <STDIN>
# SN, Summary numbers:
#   number of records   .. number of data rows in the VCF
#   number of no-ALTs   .. reference-only sites, ALT is either "." or identical to REF
#   number of SNPs      .. number of rows with a SNP
#   number of MNPs      .. number of rows with a MNP, such as CC>TT
#   number of indels    .. number of rows with an indel
#   number of others    .. number of rows with other type, for example a symbolic allele or
#                          a complex substitution, such as ACT>TCGA
#   number of multiallelic sites     .. number of rows with multiple alternate alleles
#   number of multiallelic SNP sites .. number of rows with multiple alternate alleles, all SNPs
# 
#   Note that rows containing multiple types will be counted multiple times, in each
#   counter. For example, a row with a SNP and an indel increments both the SNP and
#   the indel counter.
# 
# SN    [2]id   [3]key  [4]value
SN      0       number of samples:      1
SN      0       number of records:      4991
SN      0       number of no-ALTs:      0
SN      0       number of SNPs: 4950
SN      0       number of MNPs: 46
SN      0       number of indels:       80
SN      0       number of others:       12
SN      0       number of multiallelic sites:   121
SN      0       number of multiallelic SNP sites:       26
# TSTV, transitions/transversions:
# TSTV  [2]id   [3]ts   [4]tv   [5]ts/tv        [6]ts (1st ALT) [7]tv (1st ALT) [8]ts/tv (1st ALT)
TSTV    0       1618    3255    0.50    1608    3249    0.49
# SiS, Singleton stats:
# SiS   [2]id   [3]allele count [4]number of SNPs       [5]number of transitions        [6]number of transversions      [7]number of indels     [8]repeat-consistent    [9]repeat-inconsistent  [10]not applicable
SiS     0       1       4859    1614    3245    0       0       0       0
```

### Sample wise stats

_This section work in progress_

```math
\begin{aligned}
Sensitivity  & = \frac{TP}{TP+FN} \\
 \\
Specificity & = \frac{TN}{TN+FP} \\
\end{aligned}
```

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
_Please explain the options_

You can find the result out pgge folder [here](https://github.com/nuzla/Pangenome-Graphs-Workshop/tree/main/pgge)

You can also run this task as a SLURM job using the script [pgge_slurm.sh](https://github.com/nuzla/Pangenome-Graphs-Workshop/blob/main/Scripts/pgge_slurm.sh)

The output graph from `pgge` is shown below. 

![pgge-l100000-s5000.pgge.tsv.png](https://github.com/nuzla/Pangenome-Graphs-Workshop/blob/main/pgge/pgge-l100000-s5000.pgge.tsv.png?raw=true)
-->
