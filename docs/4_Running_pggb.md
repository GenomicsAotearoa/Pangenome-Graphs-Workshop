# 4. Running PGGB
An indtroduction to pangenome graph construction using PGGB.
If you get stuck on what a concept is or command does please refer to the extensive [documentation](https://pggb.readthedocs.io/en/latest/index.html).


## Learning objectives
!!! info "Objectives"
    - build pangenome graphs using pggb
    - explore pggb’s results
    - understand how parameters affect the built pangenome graphs
    

## How does graph building work in PGGB?

!!! info "PGGB 101"

    PBBG is a sequence graph and builds pangenome variation graphs from a set of input sequences.
    
    A pangenome variation graph is a kind of generic multiple sequence alignment. It lets us understand any kind of sequence variation between a collection of genomes. 
    Regions where genomes are similar are reflected as the same PATH through the graph, NODES represent the differences/variants between genomes.
    
    PGGB generates this kind of graph using a series of steps that each use a particular tool. 
    1. an all-to-all alignment of input sequences (wfmash)
    2. graph induction (seqwish)
    3. progressive normalization (smoothxg, gfaffix)
    
    By setting parameters PGGB can genrate several other outputs after construction of the graph.
    -m PGGB generates diagnostic visualizations of the graph (odgi)
    -V A variant call report (in VCF) representing both small and large variants in the variant graph (vg)
    
    PGGB writes its output in GFAv1 format, which can be used as input by numerous "genome graph" and pangenome tools, such as the vg and odgi toolkits.
    
    PGGB has been tested at scale in the Human Pangenome Reference Consortium (HPRC) as a method to build a graph from the draft human pangenome. 
    More details can be find [(PGGB)](https://github.com/pangenome/pggb)
    


## Getting started
NeSI HPC environment is used in this workshop for the analysis but you should be able to reproduce this environment by following the step in the (Tools and Setup)[2_Tools_and_setup.md] section. 
The following sections assume that you have created the correct folder structure and copied over the approriate data.

### Preparing and checking input sequences

Create an index for the sequence using SAMtools and check.

!!! terminal "code"

    ```bash    
    module purge
    module load SAMtools/1.16.1-GCC-11.3.0
    samtools faidx XXX.fa
    ```

    Inspect the index.

    ```bash
    tail XXX.fa.fai
    ```

    !!! success "Output"
        ```
        XXX       2248966 16      60      61
        YYY       2249014 2286474 2249014 2249015
        ZZZ       2248965 4535499 2248965 2248966
        QQQ       2249050 6784474 2249050 2249051
        ```

### Testing the PGGB tool

!!! terminal "code"

    ```bash
    module purge
    module load pggb/0.5.3-Miniconda3

    # Execute `pggb --help` to check the command list of PGGB.
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

### Key parameters for executing [PGGB](https://github.com/pangenome/pggb)

The overall structure of pggb's output graph is defined by three parameters: genome number (-n), segment length (-s), and pairwise identity (-p). 

Genome number (-n) is a given, but varies in ways that are difficult to infer and is thus left up to the user. Segment length defines the seed length used by the "MashMap3" homology mapper in wfmash. 

The pairwise identity (-p) is the minimum allowed pairwise identity between seeds, which is estimated using a mash-type approximation based on k-mer Jaccard. Mappings are initiated from collinear chains of around 5 seeds (-l, --block-length), and extended greedily as far as possible, allowing up to -n minus 1 mappings at each query position.

An additional parameter, -k, can also greatly affect graph structure by pruning matches shorter than a given threshold from the initial graph model. In effect, -k N removes any match shorter than Nbp from the initial alignment. This filter removes potentially ambiguous pairwise alignments from consideration in establishing the initial scaffold of the graph.

The initial graph is defined by parameters to wfmash and seqwish. But due to the ambiguities generated across the many pairwise alignments we use as input, this graph can be locally very complex. To regularize it we orchestrate a series of graph transformations. First, with smoothxg, we "smooth" it by locally realigning sequences to each other with a traditional multiple sequence alignment (we specifically apply POA). This process repeats multiple times to smooth over any boundary effects that may occur due to binning errors near MSA boundaries. Finally, we apply gfaffix to remove forks where both alternatives have the same sequence.

-S generate statistics of the seqwish and smoothxg graph

-m generate MultiQC report of graphs' statistics and visualizations, automatically runs odgi stats

-V specify a set of VCFs to produce with SPEC = REF:DELIM[:LEN][,REF:DELIM:[LEN]]* the paths matching ^REF are used as a reference, while the sample haplotype are derived from path names, e.g. when DELIM=# and with '-V chm13:#', a path named HG002#1#ctg would be assigned to sample HG002 phase 1. If LEN is specified and greater than 0, the VCFs are decomposed, filtering sites whose max allele length is greater than LEN. [default: off]

-o, --output-dir PATH       output directory

### Examples of key parameters for executing PGGB
- Human, whole genome, 90 haplotypes:
  - ```bash pggb -p 98 -s 50k -n 90 -k 79 ...  ```
- 15 helicobacter genomes, 5% divergence:
  - ```bash pggb -n 15 -k 79, and 15 at higher (10%) divergence pggb -n 15 -k 19 -P asm20 ...  ```
- Yeast genomes, 5% divergence: pggb's defaults should work well, just set -n.
- Aligning 9 MHC class II assemblies from vertebrate genomes (5-10% divergence):
  - ```bash pggb -n 9 -k 29 ... ```
- A few thousand bacterial genomes
  - ```bash pggb -x auto -n 2146 .... In general mapping sparsification (-x auto) is a good idea when you have many hundreds to thousands of genomes. ```
- pggb defaults to using the number of threads as logical processors on the system (the thread count given by getconf _NPROCESSORS_ONLN).
  - Use -t to set an appropriate level of parallelism if you can't use all the processors on your system.

### Use mash triangle to check the pairwise identity of the input genomes, which will give us some idea how to set -p 

!!! terminal "code"

    ```bash
    module purge
    module load Mash/2.3-GCC-11.3.0

    mash triangle XXX.fa > XXX.fa_mash
    
    # Inspect the output
    tail XXX.fa_mash
    ```

    ??? success "Output"

        ```
                5
        XXX
        YYY    0.0010072
        ZZZ    0.00121124      0.000830728
        QQQ    0.00251903      0.00366686      0.00375609
        RRR    
        ```


<b> NEEDS EXPLANATION OF OUTPUT </b>

## Running PGGB

### Generate a pangenome graph
We will be usign the 5 genomes to generate a pangenomegraph.




### Construct pangenome graph for 4Sim genomes with -k 1000, -p 96

!!! terminal "code"

     ```bash
     module purge
     module load pggb/0.5.3-Miniconda3
     
     # Execute pggb, set -s 1000
     pggb -i 4Sim.fa -s 1000 -p 96 -n 4 -t 24 -S -m -o 4Sim_1K96 -V 'NC_017518:#'
     ```

<hr>
### Extened learning: Running `pggb` as a [SLURM](https://github.com/SchedMD/slurm) Job

<b> Please do NOT run the code below, this is an example for power users </b>

Executing shell scripts in the NeSI environment might not be the best way to handle larger files which will require large memory, CPU power and time. 
We can modify the previously explained script as below to run as SLURM job. Note the additional parameters specified by `#SBATCH` which will indicate maximum resource limitations. 

The following is a SLURM script (`pggb_4Sim.sl`) for 3 PGGB runs with different settings

    ```bash
    #!/bin/bash -e     
    #SBATCH --account       nesi02659
    #SBATCH --job-name      pggb_4Sim
    #SBATCH --cpus-per-task 16
    #SBATCH --mem           16G
    #SBATCH --time          1:00:00
    #SBATCH --error         %x_%j.err
    #SBATCH --output        %x_%j.out
    
    # Modules
    module purge
    module load pggb/0.5.3-Miniconda3    
    
    # Variables
    WD=~/pg_workshop #Working Directory
    data=${WD}/4Sim.fa  
    
    # Run PGGB
    # 1K96
    pggb -i $data -s 1000 -p 96 -n 4 -t $SLURM_CPUS_PER_TASK -S -m -o $WD/4Sim_1K96 -V 'NC_017518:#'
    # 10K96
    pggb -i $data -s 10000 -p 96 -n 4 -t $SLURM_CPUS_PER_TASK -S -m -o $WD/4Sim_10K96 -V 'NC_017518:#'
    # 10K96_K79
    pggb -i $data -s 1000 -p 96 -n 4 -K 79 -t $SLURM_CPUS_PER_TASK -S -m -o $WD/4Sim_1K96_K79 -V 'NC_017518:#'
    ```

The job can be submitted using the `sbatch` command as follows. Take note of the job ID for tracking the run.

    ```bash
    sbatch pggb_4Sim.sl
    ```
