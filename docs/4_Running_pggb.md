# 4. Running PGGB

### Learning objectives

!!! quote ""

    Build pangenome graphs using PGGB


## Getting started
NeSI HPC environment is used for the analysis. Please make sure you have a NeSI account and you are able to login.

## Construct pangenome graph for the five *Neisseria meningitidis* genomes

Create an index for the sequence using SAMtools and check.

!!! terminal "code"

    ```bash    
    module purge
    module load SAMtools/1.16.1-GCC-11.3.0
    samtools faidx 5NM.fa
    ```

Inspect the index.

!!! terminal "code"

    ```bash
    more 5NM.fa.fai
    ```

    !!! success "Output"
        
        ```
        NC_003112.2     2272360   60        80      81
        NC_017518.1     2248966   2300889   80      81
        NZ_CP007668.1   2324822   4578040   80      81
        NZ_CP016880.1   2207174   6931986   80      81
        NZ_CP020423.2   2244886   9166836   80      81
        ```

### Executing `pggb` 

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



## Running PGGB

### Use `mash triangle` to check the pairwise identity of the input genomes, which will give us some idea how to set `-p` 

!!! terminal "code"

    ```bash
    module purge
    module load Mash/2.3-GCC-11.3.0

    mash triangle 5NM.fa > 5NM.fa_mash
    
    # Inspect the output
    more 5NM.fa_mash
    ```

    ??? success "Output"

        ```
                5
        NC_003112.2
        NC_017518.1     0.0152404
        NZ_CP007668.1   0.0149234       0.00635099
        NZ_CP016880.1   0.0178909       0.0171265       0.0170111
        NZ_CP020423.2   0.0190552       0.0194352       0.0185579       0.0106974
        ```


<b>The lower triangle represents the pairwise distances between the 5NM genomes. We can observe that the largest paired distance is 0.0190552, which is approximately 0.02. Considering that lower values indicate better alignment, we are going to use an alignment threshold of `-p 94` for constructing the pangenome graph. </b>

### Construct pangenome graph for 5NM genomes with `-k 1000`, `-p 96`

!!! terminal "code"

     ```bash
     module purge
     module load pggb/0.5.3-Miniconda3
     
     # Execute pggb, set -s 2000 and -p 94
     pggb -i 5NM.fa -s 2000 -p 94 -n 4 -t 24 -S -m -o 5NM_2K94 -V 'NC_017518.1:#'
     ```

<hr>
### Extened learning: Running `pggb` as a [SLURM](https://github.com/SchedMD/slurm) Job

<b> Please do NOT run the code below, this is an example for power users </b>

Executing shell scripts in the NeSI environment might not be the best way to handle larger files which will require large memory, CPU power and time. 
We can modify the previously explained script as below to run as SLURM job. Note the additional parameters specified by `#SBATCH` which will indicate maximum resource limitations. 

The following is a SLURM script (`pggb_5NM_2k94.sl`) for PGGB with -k 2000, and -p 94 

    ```bash
    #!/bin/bash -e     
    #SBATCH --account       nesi02659
    #SBATCH --job-name      pggb_5NM
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
    data=${WD}/5NM.fa  
    
    # Run PGGB
    # 2K94
    pggb -i $data -s 2000 -p 94 -n 4 -t $SLURM_CPUS_PER_TASK -S -m -o $WD/5NM_2K94 -V 'NC_017518.1:#'
    
    ```

The job can be submitted using the `sbatch` command as follows. Take a note of the job ID for tracking the run.

    ```bash
    sbatch pggb_5NM.sl
    ```
