# pangenome graph construction with PGGB

In this workshop, we employes PanGenome Graph Builder (PGGB) to construction the pangenome graphs

## How does the pggb graph build work?

!!! info ""

    pggb builds pangenome variation graphs from a set of input sequences.
    
    A pangenome variation graph is a kind of generic multiple sequence alignment. It lets us understand any kind of sequence variation between a collection of genomes. It shows us similarity where genomes walk through the same parts of the graph, and differences where they do not.
    
    pggb generates this kind of graph using an all-to-all alignment of input sequences (wfmash), graph induction (seqwish), and progressive normalization (smoothxg, gfaffix). After construction, pggb generates diagnostic visualizations of the graph (odgi). A variant call report (in VCF) representing both small and large variants can be generated based on any reference genome included in the graph (vg). pggb writes its output in GFAv1 format, which can be used as input by numerous "genome graph" and pangenome tools, such as the vg and odgi toolkits.
    
    pggb has been tested at scale in the Human Pangenome Reference Consortium (HPRC) as a method to build a graph from the draft human pangenome. 
    
    more details can be find [(PGGB)](https://github.com/pangenome/pggb)
    
### Learning objectives

!!! quote ""

    - build pangenome graphs using pggb
    - explore pggb’s results
    - understand how parameters affect the built pangenome graphs

## Getting started
NeSI HPC environment is used for the analysis. Please make sure to have a NeSI account and you are able to login.

### Setting up your project directory and download the datasets

!!! terminal "code"

    ```bash
    # Create a new directory in somewhere and change to that directory
    mkdir pg_workshop
    cd pg_workshop
    # Keep a note of the absolute path of your directory
    pwd
    /home/zyang/pg_worhshop
    
    # Downloading and preparing datasets
    git clone https://github.com/ZoeYang2020/dataset_for_pg_workshop
    
    # copy the 4Sim.fa dataset to your work directory, mine is /home/zyang/pg_worhshop
    cp /home/zyang/pg_worhshop/dataset_for_pg_workshop/datasets_for_PangenomeGraphConstruction_pg_workshop/4Sim.fa /home/zyang/pg_worhshop
    
    # go back to your work directory 
    cd /home/zyang/pg_worhshop
    ```

## Construct pangenome graph for the 4Sim genomes

!!! terminal "code"

    ```bash
    #Creating an index for the seqence file using samtools and check
    #In Nesi environment you will have to load the samtools module first
    
    module purge
    module load SAMtools/1.16.1-GCC-11.3.0
    samtools faidx 4Sim.fa
    less -S 4Sim.fa.fai
    NC_017518       2248966 16      60      61
    ST41Sim 2249014 2286474 2249014 2249015
    ST154Sim        2248965 4535499 2248965 2248966
    ST42Sim 2249050 6784474 2249050 2249051
    ```
### Executing `pggb` tool

!!! terminal "code"

    ```bash
    module purge
    module load pggb/0.5.3-Miniconda3
    

    # Execute `pggb --help` to check the command list of PGGB
    pggb --help
    ```
    ??? success "Output"
        ```bash 
        ERROR: mandatory arguments -i and -n
        ERROR: -n must be greater than or equal to 2
        usage: /usr/local/bin/pggb -i <input-fasta> -n <n-haplotypes> [options]
        options:
           [wfmash]
            -i, --input-fasta FILE      input FASTA/FASTQ file
            -s, --segment-length N      segment length for mapping [default: 5000]
            -l, --block-length N        minimum block length filter for mapping [default: 5*segment-length]
            -p, --map-pct-id PCT        percent identity for mapping/alignment [default: 90]
            -n, --n-haplotypes N        number of haplotypes
            -N, --no-split              disable splitting of input sequences during mapping [default: enabled]
            -x, --sparse-map N          keep this fraction of mappings ('auto' for giant component heuristic) [default: 1.0]
            -K, --mash-kmer N           kmer size for mapping [default: 19]
            -F, --mash-kmer-thres N     ignore the top % most-frequent kmers [default: 0.001]
            -Y, --exclude-delim C       skip mappings between sequences with the same name prefix before
                                        the given delimiter character [default: all-vs-all and !self]
           [seqwish]
            -k, --min-match-len N       filter exact matches below this length [default: 19]
            -f, --sparse-factor N       keep this randomly selected fraction of input matches [default: no sparsification]
            -B, --transclose-batch      number of bp to use for transitive closure batch [default: 10000000]
           [smoothxg]
            -X, --skip-normalization    do not normalize the final graph [default: normalize the graph]
            -H, --n-haplotypes-smooth N number of haplotypes, if different than that set with -n [default: -n]
            -j, --path-jump-max         maximum path jump to include in block [default: 0]
            -e, --edge-jump-max N       maximum edge jump before breaking [default: 0]
            -G, --poa-length-target N,M target sequence length for POA, one per pass [default: 700,900,1100]
            -P, --poa-params PARAMS     score parameters for POA in the form of match,mismatch,gap1,ext1,gap2,ext2
                                        may also be given as presets: asm5, asm10, asm15, asm20
                                        [default: 1,19,39,3,81,1 = asm5]
            -O, --poa-padding N         pad each end of each sequence in POA with N*(mean_seq_len) bp [default: 0.001]
            -d, --pad-max-depth N       depth/haplotype at which we don't pad the POA problem [default: 100]
            -b, --run-abpoa             run abPOA [default: SPOA]
            -z, --global-poa            run the POA in global mode [default: local mode]
            -M, --write-maf             write MAF output representing merged POA blocks [default: off]
            -Q, --consensus-prefix P    use this prefix for consensus path names [default: Consensus_]
           [odgi]
            -v, --skip-viz              don't render visualizations of the graph in 1D and 2D [default: make them]
            -S, --stats                 generate statistics of the seqwish and smoothxg graph [default: off]
           [vg]
            -V, --vcf-spec SPEC         specify a set of VCFs to produce with SPEC = REF:DELIM[:LEN][,REF:DELIM:[LEN]]*
                                        the paths matching ^REF are used as a reference, while the sample haplotypes
                                        are derived from path names, e.g. when DELIM=# and with '-V chm13:#',
                                        a path named HG002#1#ctg would be assigned to sample HG002 phase 1.
                                        If LEN is specified and greater than 0, the VCFs are decomposed, filtering
                                        sites whose max allele length is greater than LEN. [default: off]
           [multiqc]
            -m, --multiqc               generate MultiQC report of graphs' statistics and visualizations,
                                        automatically runs odgi stats [default: off]
           [general]
            -o, --output-dir PATH       output directory
            -D, --temp-dir PATH         directory for temporary files
            -a, --input-paf FILE        input PAF file; the wfmash alignment step is skipped
            -r, --resume                do not overwrite existing outputs in the given directory
                                        [default: start pipeline from scratch]
            -t, --threads N             number of compute threads to use in parallel steps [default: 72]
            -T, --poa-threads N         number of compute threads to use during POA (set lower if you OOM during smoothing)
            -A, --keep-temp-files       keep intermediate graphs
            -Z, --compress              compress alignment (.paf), graph (.gfa, .og), and MSA (.maf) outputs with pigz,
                                        and variant (.vcf) outputs with bgzip
            --version                   display the version of pggb
            -h, --help                  this text
        
        Use wfmash, seqwish, smoothxg, odgi, gfaffix, and vg to build, project and display a pangenome graph.
        ```
### key parameters for executing PGGB
https://github.com/pangenome/pggb
The overall structure of pggb's output graph is defined by three parameters: genome number (-n), segment length (-s), and pairwise identity (-p). 

Genome number (-n) is a given, but varies in ways that are difficult to infer and is thus left up to the user. Segment length defines the seed length used by the "MashMap3" homology mapper in wfmash. 

The pairwise identity (-p) is the minimum allowed pairwise identity between seeds, which is estimated using a mash-type approximation based on k-mer Jaccard. Mappings are initiated from collinear chains of around 5 seeds (-l, --block-length), and extended greedily as far as possible, allowing up to -n minus 1 mappings at each query position.

An additional parameter, -k, can also greatly affect graph structure by pruning matches shorter than a given threshold from the initial graph model. In effect, -k N removes any match shorter than Nbp from the initial alignment. This filter removes potentially ambiguous pairwise alignments from consideration in establishing the initial scaffold of the graph.

The initial graph is defined by parameters to wfmash and seqwish. But due to the ambiguities generated across the many pairwise alignments we use as input, this graph can be locally very complex. To regularize it we orchestrate a series of graph transformations. First, with smoothxg, we "smooth" it by locally realigning sequences to each other with a traditional multiple sequence alignment (we specifically apply POA). This process repeats multiple times to smooth over any boundary effects that may occur due to binning errors near MSA boundaries. Finally, we apply gfaffix to remove forks where both alternatives have the same sequence.



### examples of key parameters for executing PGGB
- Human, whole genome, 90 haplotypes: pggb -p 98 -s 50k -n 90 -k 79 ...
- 15 helicobacter genomes, 5% divergence: pggb -n 15 -k 79, and 15 at higher (10%) divergence pggb -n 15 -k 19 -P asm20 ...
- Yeast genomes, 5% divergence: pggb's defaults should work well, just set -n.
- Aligning 9 MHC class II assemblies from vertebrate genomes (5-10% divergence): pggb -n 9 -k 29 ...
- A few thousand bacterial genomes pggb -x auto -n 2146 .... In general mapping sparsification (-x auto) is a good idea when you have many hundreds to thousands of genomes.
- pggb defaults to using the number of threads as logical processors on the system (the thread count given by getconf _NPROCESSORS_ONLN). Use -t to set an appropriate level of parallelism if you can't use all the processors on your system.


### other parameters for executing PGGB
-S generate statistics of the seqwish and smoothxg graph

-m generate MultiQC report of graphs' statistics and visualizations, automatically runs odgi stats

-V specify a set of VCFs to produce with SPEC = REF:DELIM[:LEN][,REF:DELIM:[LEN]]* the paths matching ^REF are used as a reference, while the sample haplotype are derived from path names, e.g. when DELIM=# and with '-V chm13:#', a path named HG002#1#ctg would be assigned to sample HG002 phase 1. If LEN is specified and greater than 0, the VCFs are decomposed, filtering sites whose max allele length is greater than LEN. [default: off]

-o, --output-dir PATH       output directory

### Use mash triangle to check the pairwise identity of the input genomes, which will give us some idea how to set -p 

!!! terminal "code"

    ```bash
    module purge
    module load Mash/2.3-GCC-11.3.0

    mash triangle 4Sim.fa >4Sim.fa_mash
    ```
    ```bash
    less -S 4Sim.fa_mash
            4
    NC_017518
    ST41Sim 0.0010072
    ST154Sim        0.00121124      0.000830728
    ST42Sim 0.00251903      0.00366686      0.00375609
    ```


### construct pangenome graph for 4Sim genomes with Singularity container PGGB, -k 1000, -p 96

!!! terminal "code"

     ```bash
     module purge
     module load pggb/0.5.3-Miniconda3
     
     # Execute singularity exec ${container} pggb, set -s 1000
     pggb -i 4Sim.fa -s 1000 -p 96 -n 4 -t 24 -S -m -o 4Sim_1K96 -V 'NC_017518:#'
     ```
### Executing `pggb` as a [SLURM](https://github.com/SchedMD/slurm) Job
Executing shell scripts in the Nesi environment might not be the best way to handle larger files which will require large memory, CPU power and time. We can modify the previusely explained script as below to run as SLURM job. Note the additional parameters specified by `#SBATCH` which will indicate maximum resource limitations. 



#### pggb_slurm_1K96_4Sim.sh

!!! terminal "code"

    ```bash
    #!/bin/bash
    
    #SBATCH --account       ga03793
    #SBATCH --job-name      4Sim_1K96
    #SBATCH --cpus-per-task 8
    #SBATCH --mem           4G
    #SBATCH --time          1:00:00
    
    module purge
    module load pggb/0.5.3-Miniconda3
    
    #export container to a variable for convenience
    WD=/home/zyang/pg_workshop #Working Directory
    data=/home/zyang/pg_workshop/4Sim.fa
    output=/home/zyang/pg_workshop
    
    
    pggb -i $data -s 1000 -p 96 -n 4 -t 24 -S -m -o $output/4Sim_1K96 -V 'NC_017518:#'   
    ```
The job can be submitted using the `sbatch` command it will show a job id.


#### pggb_slurm_10K96_4Sim.sh

!!! terminal "code"

    ```bash
    #!/bin/bash
    
    #SBATCH --account       ga03793
    #SBATCH --job-name      4Sim_10K96
    #SBATCH --cpus-per-task 8
    #SBATCH --mem           4G
    #SBATCH --time          1:00:00
    
    module purge
    module load pggb/0.5.3-Miniconda3
    
    #export container to a variable for convenience
    WD=/home/zyang/pg_workshop #Working Directory
    data=/home/zyang/pg_workshop/4Sim.fa
    output=/home/zyang/pg_workshop

    pggb -i $data -s 10000 -p 96 -n 4 -t 24 -S -m -o $output/4Sim_10K96 -V 'NC_017518:#'
    ```
#### pggb_slurm_10K96_K79_4Sim.sh

!!! terminal "code"

    ```bash
    #!/bin/bash
    
    #SBATCH --account       ga03793
    #SBATCH --job-name      4Sim_1K96_K79
    #SBATCH --cpus-per-task 8
    #SBATCH --mem           4G
    #SBATCH --time          1:00:00
    
    module purge
    module load Singularity
    
    #export container to a variable for convenience
    WD=/home/zyang/pg_workshop #Working Directory
    container=/nesi/project/nesi02659/software/pggb/pggb_0.5.3.simg
    data=/home/zyang/pg_workshop/4Sim.fa
    output=/home/zyang/pg_workshop
    
    
    #Bind filesystem to container image
    export SINGULARITY_BIND="${WD}, /nesi/project/nesi02659/"
    
    singularity exec ${container} pggb -i $data -s 1000 -p 96 -n 4 -K 79 -t 24 -S -m -o $output/4Sim_1K96_K79 -V 'NC_017518:#'
    ```

### Evaluate Pangenome Graphs for 4Sim Genomes Constructed with Different Settings
- We have employed three distinct settings to construct the pangenome graph of the 4Sim genomes. Which setting yielded the most optimal result? How can we determine this? 

- download the multiqc.html file, check the Detailed ODGI stats table.
#### 1k96
| Sample Name                         | Length    | Nodes  | Edges  | Components | A   |C    |T    |G    |N   |
|:-----                               |----------:|-------:|-------:|-----------:|----:|----:|----:|----:|----:|
|seqwish	|2280344	|22216	|29823	|4	|1	|551639	|578590	|557450	|592665	|0|
|smooth	|2261163	|29965	|40179	|4	|1	|548754	|574693	|551650	|586066	|0|

![1k96 ODGI 1D visualization by path orientation](https://github.com/ZoeYang2020/Pangenome-Graphs-Workshop/blob/main/pggb_4sim/1k96/4Sim.fa.97e7156.417fcdf.7659dc8.smooth.final.og.viz_inv_multiqc.png?raw=true])


#### 1k96,-K79


| Sample Name                         | Length    | Nodes  | Edges  | Components | A   |C    |T    |G    |N   |
|:-----                               |----------:|-------:|-------:|-----------:|----:|----:|----:|----:|----:|
|seqwish	|2279905	|22209	|29812	|4	|1	|551588	|578459	|557291	|592567	|0|
|smooth	|2261401	|29976	|40199	|4	|1	|553049	|580469	|547460	|580423	|0|


![1k96K79 ODGI 1D visualization by path orientation](https://github.com/ZoeYang2020/Pangenome-Graphs-Workshop/blob/main/pggb_4sim/1k96_K79/4Sim.fa.f958389.417fcdf.7659dc8.smooth.final.og.viz_inv_multiqc.png?raw=true])

#### 10k96
| Sample Name                         | Length    | Nodes  | Edges  | Components | A   |C    |T    |G    |N   |
|:-----                               |----------:|-------:|-------:|-----------:|----:|----:|----:|----:|----:|
|seqwish	|2340700	|22166	|29759	|4	|1	|566559	|594741	|571836	|607564	|0|
|smooth	|2319601	|29888	|40070	|4	|1	|566755	|599287	|562078	|591481	|0|

![10k96 ODGI 1D visualization by path orientation](https://github.com/ZoeYang2020/Pangenome-Graphs-Workshop/blob/main/pggb_4sim/10k96/4Sim.fa.e7f7fe6.417fcdf.7659dc8.smooth.final.og.viz_inv_multiqc.png?raw=true])


### vg deconstruct graph to get the variations in vcf 

```bash
mkdir vg_deconstruct
#copy gfa to vg_deconsturct and rename the gfa files with its PGGB settings
cp /home/zyang/pg_workshop/4Sim_1K96/4Sim.fa.97e7156.417fcdf.7659dc8.smooth.final.gfa /home/zyang/pg_workshop/vg_deconstruct/4Sim_1K96.gfa
cp /home/zyang/pg_workshop/4Sim_1K96_K79/4Sim.fa.f958389.417fcdf.7659dc8.smooth.final.gfa /home/zyang/pg_workshop/vg_deconstruct/4Sim_1K96_K79.gfa

cd /home/zyang/pg_workshop/vg_deconstruct

module purge
module load vg/1.46.0
module load BCFtools/1.15.1-GCC-11.3.0

vg deconstruct -p NC_017518  -a -e 4Sim_1K96.gfa >4Sim_1K96_aep1.vcf
bcftools stats 4Sim_1K96_aep1.vcf >4Sim_1K96_aep1.vcf_stats
```

#### 4Sim_vg_deconstruct.sh
```bash
#!/usr//bin/bash

#SBATCH --account       ga03793
#SBATCH --job-name      4Sim_vg_deconstruct
#SBATCH --cpus-per-task 8
#SBATCH --mem           4G
#SBATCH --time          1:00:00

module purge
module load vg/1.46.0
module load BCFtools/1.15.1-GCC-11.3.0

export container to a variable for convenience
inputGFA=/home/zyang/pg_workshop/vg_deconstruct/*.gfa
input_folder=/home/zyang/pg_workshop/vg_deconstruct
output=/home/zyang/pg_workshop/vg_deconstruct

for f in $inputGFA
do

x=$(basename $f .gfa)
echo ${x}

vg deconstruct -p NC_017518  -a -e $input_folder/${x}.gfa > $output/${x}aep1.vcf
bcftools stats $output/${x}aep1.vcf >$output/${x}aep1.vcf_stats
done
```
#### bcftools isec to check the overlap of the 1k96 -K 19, 1k96, -K 79
```bash
#!/bin/bash

#SBATCH --account       ga03793
#SBATCH --job-name      4Sim_1k96_isec
#SBATCH --cpus-per-task 8
#SBATCH --mem           4G
#SBATCH --time          1:00:00

module purge
module load BCFtools/1.15.1-GCC-11.3.0

input_folder=/home/zyang/pg_workshop/vg_deconstruct
output_folder=/home/zyang/pg_workshop/vg_deconstruct


bcftools view $input_folder/4Sim_1K96aep1.vcf  -Oz -o $output_folder/4Sim_1K96aep1.vcf.gz
bcftools view $input_folder/4Sim_1K96_K79aep1.vcf -Oz -o $output_folder/4Sim_1K96_K79aep1.vcf.gz

bcftools index $output_folder/4Sim_1K96aep1.vcf.gz
bcftools index $output_folder/4Sim_1K96_K79aep1.vcf.gz

bcftools isec $output_folder/4Sim_1K96aep1.vcf.gz $output_folder/4Sim_1K96_K79aep1.vcf.gz -p $output_folder/isec_4Sim_1K96
```

```bash
#specific to  1k96 -K 19
less -S /home/zyang/pg_workshop/vg_deconstruct/isec_4Sim_1K96/0000.vcf
#specific to  1k96 -K 79
less -S /home/zyang/pg_workshop/vg_deconstruct/isec_4Sim_1K96/0001.vcf
```
##### difference between 1k96 -K 19 Vs 1k96 -K 79
![difference between 1k96 -K 19 Vs 1k96 -K 79](https://github.com/ZoeYang2020/Pangenome-Graphs-Workshop/blob/main/Figures/4Sim1K96_K19vsK79_specific_variation.png?raw=true)


### PGGE 
This pangenome graph evaluation pipeline measures the reconstruction accuracy of a pangenome graph (in the variation graph model). Its goal is to give guidance in finding the best pangenome graph construction tool for a given input data and task.

```bash
mkdir 4Sim_pgge
#copy *.gfa to 4Sim_pgge
cp /home/zyang/pg_workshop/vg_deconstruct/*.gfa /home/zyang/pg_workshop/4Sim_pgge
cp /home/zyang/pg_workshop/4Sim_10K96/4Sim.fa.e7f7fe6.417fcdf.7659dc8.smooth.final.gfa /home/zyang/pg_workshop/4Sim_pgge/4Sim_10K96.gfa
```
#### PGGE script
```bash
#!/bin/bash

#SBATCH --account       ga03793
#SBATCH --job-name      4Sim_pgge
#SBATCH --cpus-per-task 8
#SBATCH --mem           4G
#SBATCH --time          1:00:00

module purge
module load Singularity

#export container to a variable for convenience

container=/nesi/project/nesi02659/software/pgge/pgge_032023.simg


WD=/home/zyang/pg_workshop/4Sim_pgge #Working Directory
inputGFA=/home/zyang/pg_workshop/4Sim_pgge/*.gfa
input_folder=/home/zyang/pg_workshop/4Sim_pgge
inputfa=/home/zyang/pg_workshop/4Sim.fa
output=/home/zyang/pg_workshop/4Sim_pgge
beehave=/home/zyang/pg_workshop/beehave.R


#Bind filesystem to container image
export SINGULARITY_BIND="${WD}, /nesi/project/nesi02659/"


for x in $inputGFA

do
singularity exec ${container} pgge -g $x -f $inputfa -o $output -r $beehave -b $output/pgge_4Sim_peanut_bed -l 100000 -s 5000 -t 16

done
```

## ODGI paths to extract distance  
```bash
mkdir odgi_distance
cp odgi_distance
#copy gfa to odgi_distance working directory 
cp /home/zyang/pg_workshop/vg_deconstruct/4Sim_1K96.gfa /home/zyang/pg_workshop/odgi_distance/

module purge
module load Singularity
container=/nesi/project/nesi02659/software/odgi/odgi_0.8.2.simg
singularity exec ${container} odgi paths -i 4Sim_1K96.gfa -d -D 'AAAA' >4Sim_1K96.gfa_distance
cut -f 1,2,6 4Sim_1K96.gfa_distance >4Sim_1K96.gfa_distance_cut
```
### script for ODGI paths to extract distance
```bash
#!/bin/bash

#SBATCH --account       ga03793
#SBATCH --job-name      4Sim_1K96_distance
#SBATCH --cpus-per-task 8
#SBATCH --mem           4G
#SBATCH --time          1:00:00

module purge
module load Singularity

#export container to a variable for convenience
container=/nesi/project/nesi02659/software/odgi/odgi_0.8.2.simg
data=/home/zyang/pg_workshop/odgi_distance/4Sim_1K96.gfa
output=/home/zyang/pg_workshop/odgi_distance

singularity exec ${container} odgi paths -i $data -d -D 'AAAA' >$output/4Sim_1K96.gfa_distance
cut -f 1,2,6 $output/4Sim_1K96.gfa_distance >$output/4Sim_1K96.gfa_distance_cut
```
### R script for clustering based on distance among paths of a graph, 4Sim 1k96

list2dist_clustering_4Sim.R
```bash
setwd("/home/zyang/pg_workshop/odgi_distance")

library(reshape)
library(ape)

# read in the data
dat=read.csv("4Sim_1K96.gfa_distance_cut",sep="\t")
dat
# use reshape's cast function to change to matrix
m <- cast(dat, group.a ~ group.b)
m
# set the row names
rownames(m) <- m[,1]
rownames(m)

#The fellowing two lines code will cause no IDs in the clustering
# get rid of a couple of rows
#m <- m[,-2]
# convert any 0s that were read in as strings to integers
#m <- apply(m, 2, as.numeric )
m

# change the matrix to a distance matrix
d <- dist(m)
d

# do hierarchical clustering
h <- hclust(d)

h
# plot the dendrogram
plot(h)

# use ape's as phylo function
tree <- as.phylo(h)
# export as newick for viewing in figtree
write.tree(phy=tree, file = '4Sim_1k96_distance.tree')
```
### run the R script for clustering based on distance on Nesi
```bash
#!/bin/bash

#SBATCH --account       ga03793
#SBATCH --job-name      4Sim_1K96_distance_clustering
#SBATCH --cpus-per-task 8
#SBATCH --mem           4G
#SBATCH --time          1:00:00

module purge
module load R/4.0.1-gimkl-2020

Rscript 8_list2dist_clustering_4Sim.R
```

## Construct pangenome graph for the 3ST genomes
### prepare dataset and build index
```bash
# copy the 3ST.fa dataset to your work directory, mine is /home/zyang/pg_worhshop
cp /home/zyang/pg_workshop/dataset_for_pg_workshop/datasets_for_PangenomeGraphConstruction_pg_workshop/3ST.fa /home/zyang/pg_workshop

# go back to your work directory 
cd /home/zyang/pg_worhshop

#build index for 3ST.fa

module purge
module load SAMtools/1.16.1-GCC-11.3.0
samtools faidx 3ST.fa

#check index 
less -S 3ST.fa.fai
NC_017518       2248966 77      60      61
ST41    2217832 2286541 60      61
ST154   2233582 4541354 60      61

```

### Use mash triangle to check the pairwise identity of the input genomes, which will give us some idea how to set -p 
```bash
module purge
module load Mash/2.3-GCC-11.3.0
mash triangle 4Sim.fa >4Sim.fa_mash
less -S 3ST.fa_mash
        3
NC_017518
ST41    0.00146992
ST154   0.00165343      0.00131423
```

### check the Mauve aligment of 3ST
Mauve alignments demonstrated large inversions among the 3ST genomes. 
![Mauve alignment of the 3STs genomes](https://github.com/ZoeYang2020/Pangenome-Graphs-Workshop/blob/main/Figures/Fig.3ST_mauve%20alignment.png??raw=true])

### pggb_slurm_2K95_3ST.sh, -k 2000, -p 95
```bash
#!/bin/bash

#SBATCH --account       ga03793
#SBATCH --job-name      3ST_2K95
#SBATCH --cpus-per-task 8
#SBATCH --mem           4G
#SBATCH --time          1:00:00

module purge
module load Singularity

#export container to a variable for convenience
container=/nesi/project/nesi02659/software/pggb/pggb_0.5.3.simg
data=/home/zyang/pg_workshop/3ST.fa
output=/home/zyang/pg_workshop

singularity exec ${container} pggb -i $data -s 2000 -p 95 -n 3 -t 24 -S -m -o $output/3ST_2K95 -V 'NC_017518:#'
```

## Construct pangenome graph for the 24NM genomes
### prepare dataset and build index
```bash
# copy the 24NM.fa.gz dataset to your work directory, mine is /home/zyang/pg_worhshop
cp /home/zyang/pg_workshop/dataset_for_pg_workshop/datasets_for_PangenomeGraphConstruction_pg_workshop/24NM.fa.gz /home/zyang/pg_workshop

# go back to your work directory 
cd /home/zyang/pg_worhshop

#uncompress 24NM.fa.gz
gzip -d 24NM.fa.gz

#build index for 24NM.fa

module purge 
module load SAMtools/1.16.1-GCC-11.3.0
samtools faidx 24NM.fa

#check index 
less -S 24NM.fa.fai
NC_003112.2     2272360 60      60      61
NC_003116.1     2184406 2310354 60      61
NC_008767.1     2194961 4531228 60      61
NC_010120.1     2153416 6762834 60      61
NC_017501.1     2277550 8952201 60      61
NC_013016.1     2145295 11267774        60      61
NC_017505.1     2242947 13448888        60      61
NC_017513.1     2184862 15729279        60      61
NC_017514.1     2223518 17950622        60      61
NC_017515.1     2250449 20211265        60      61
NC_017512.1     2227255 22499286        60      61
NZ_CP012392.1   2170619 24763743        60      61
NZ_CP016627.1   2173408 26970619        60      61
NZ_CP016646.1   2173686 29180331        60      61
NZ_CP016682.1   2175832 31390326        60      61
NZ_CP020402.2   2305818 33602508        60      61
NZ_CP031334.1   2314390 35946837        60      61
NZ_CP031324.1   2291778 38299881        60      61
NZ_CP031328.1   2223855 40629936        60      61
NZ_CP031333.1   2280611 42890936        60      61
NZ_CP021517.1   2167947 45209638        60      61
NC_017518       2248966 47413790        60      61
ST41    2217832 49700245        60      61
ST154   2233582 51955048        60      61

```

### Use mash triangle to check the pairwise identity of the input genomes, which will give us some idea how to set -p 
```bash
module purge
module load Mash/2.3-GCC-11.3.0
mash triangle 24NM.fa >24NM.fa_mash
less -S 24NM.fa_mash

        24
NC_003112.2
NC_003116.1     0.0188675
NC_008767.1     0.0174175       0.0171844
NC_010120.1     0.0184965       0.0166683       0.0166117
NC_017501.1     0.0181313       0.0181313       0.016782        0.0190552
NC_013016.1     0.0218499       0.0190552       0.0171844       0.0191812       0.0188053
NC_017505.1     0.015943        0.0191812       0.0168963       0.0187432       0.0175939       0.0212923
NC_017513.1     0.0186195       0.0175939       0.00861543      0.0167251       0.016782        0.0174175       0.0169536
NC_017514.1     0.0154006       0.0185579       0.0167251       0.0185579       0.0172424       0.0204123       0.00977265      0.0175939
NC_017515.1     0.0115313       0.0177716       0.0157245       0.0170111       0.0168391       0.0202137       0.016782        0.0155619       0.0158335
NC_017512.1     0.0192445       0.00594242      0.0177716       0.0173006       0.0184965       0.0200822       0.0202137       0.0182524       0.0199514       0.0179508
NZ_CP012392.1   0.0183741       0.0171844       0.017359        0.0172424       0.0183132       0.0204789       0.0180108       0.0170111       0.0172424       0.0162193       0.0177716
NZ_CP016627.1   0.0181313       0.0174175       0.00279888      0.0166117       0.0172424       0.0179508       0.0175349       0.00877107      0.0171265       0.0166683       0.017359
NZ_CP016646.1   0.0178909       0.0172424       0.0163866       0.0171265       0.0174761       0.0184352       0.0189925       0.0163307       0.017653        0.0165552       0.0177716
NZ_CP016682.1   0.0180108       0.0170111       0.0162193       0.0169536       0.0171844       0.0181313       0.0189299       0.0164427       0.0175939       0.0162749       0.0175349
NZ_CP020402.2   0.0196919       0.0188675       0.0177122       0.0180709       0.0134525       0.0185579       0.0171844       0.0171265       0.0170111       0.0182524       0.0199514
NZ_CP031334.1   0.0170111       0.0180709       0.0156702       0.0163866       0.0147145       0.0183132       0.0165552       0.0149759       0.0164427       0.0152936       0.0191181
NZ_CP031324.1   0.0177716       0.0189925       0.0166683       0.0181313       0.0146626       0.0193079       0.0178909       0.0161637       0.0166117       0.0160531       0.0194991
NZ_CP031328.1   0.0181313       0.0180108       0.0168391       0.0174175       0.0162749       0.0196275       0.0182524       0.0166683       0.0157789       0.0168391       0.0184965
NZ_CP031333.1   0.0117112       0.0172424       0.0158335       0.0170687       0.0168391       0.0201478       0.016782        0.0156702       0.0158335       0.00186568      0.0171844
NZ_CP021517.1   0.0178312       0.0163866       0.0171265       0.0165552       0.0180709       0.0194352       0.0178909       0.0165552       0.0174761       0.0158335       0.0170687
NC_017518       0.0152404       0.0186195       0.0166117       0.0183741       0.016782        0.0200167       0.00916565      0.0174761       0.001839        0.015508        0.0194991
ST41    0.0150813       0.0184352       0.0164427       0.0183132       0.0170687       0.0201478       0.00912584      0.017653        0.00128842      0.0158335       0.0194991       0.
ST154   0.0151872       0.0182524       0.0164989       0.0182524       0.0171265       0.0200167       0.00952763      0.0173006       0.00175922      0.0156702       0.0193714       0.

```

### pggb_slurm_2K95_3ST.sh, -k 2000, -p 95
```bash
#!/bin/bash

#SBATCH --account       ga03793
#SBATCH --job-name      24NM_10k95
#SBATCH --cpus-per-task 8
#SBATCH --mem           4G
#SBATCH --time          3:00:00

module purge
module load Singularity

#export container to a variable for convenience
WD=/home/zyang/pg_workshop #Working Directory
container=/nesi/project/nesi02659/software/pggb/pggb_0.5.3.simg
data=/home/zyang/pg_workshop/24NM.fa
output=/home/zyang/pg_workshop


#Bind filesystem to container image
export SINGULARITY_BIND="${WD}, /nesi/project/nesi02659/"

singularity exec ${container} pggb -i $data -s 10000 -p 95 -n 24 -t 24 -S -m -o $output/24NM_10K95 -V 'NC_017518:#'
```

### pggb_slurm_2K95_3ST.sh, -k 2000, -p 95, -x 
```bash
#!/bin/bash

#SBATCH --account       ga03793
#SBATCH --job-name      24NM_10k95_X
#SBATCH --cpus-per-task 8
#SBATCH --mem           4G
#SBATCH --time          3:00:00

module purge
module load Singularity

#export container to a variable for convenience
WD=/home/zyang/pg_workshop #Working Directory
container=/nesi/project/nesi02659/software/pggb/pggb_0.5.3.simg
data=/home/zyang/pg_workshop/24NM.fa
output=/home/zyang/pg_workshop


#Bind filesystem to container image
export SINGULARITY_BIND="${WD}, /nesi/project/nesi02659/"

singularity exec ${container} pggb -i $data -s 10000 -p 95 -n 24 -x auto -t 24 -S -m -o $output/24NM_10K95x -V 'NC_017518:#'
```
