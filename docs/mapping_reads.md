# Mapping Reads using `bwa mem` (Linear Method)
_Note: folder : /nesi/nobackup/nesi02659/pg_workshop/vc_exact_compare/_
We can use below script to map the reads to the reference sequnce using linear method `bwa mem`.

!!! terminal "code"

    ```bash
    #Load required modules with specific versions
    module purge
    module load BCFtools/1.15.1-GCC-11.3.0
    module load SAMtools/1.15.1-GCC-11.3.0
    module load BWA/0.7.17-GCC-11.3.0
    module load wgsim/20111017-GCC-11.3.0
    
    mkdir vc_exact_compare
    cd vc_exact_compare
    
    #Copy the below files into the folder
    cp GCF_000191525.1_ASM19152v1_genomic.fna vc_exact_compare/
    cp Simulation*.simseq.genome.fa vc_exact_compare/
    
    ls -1trhs
    ```

!!! success "Output"

    ```
    total 14M
    2.3M GCF_000191525.1_ASM19152v1_genomic.fna
    2.3M Simulation_INDEL_5000.simseq.genome.fa
    2.5M Simulation_SNP_4000_INDEL_4000_CNV_4.simseq.genome.fa
    2.3M Simulation_SNP_4000_INDEL_4000_INV_4.simseq.genome.fa
    2.3M Simulation_SNP_4000_INDEL_4000.simseq.genome.fa
    2.3M Simulation_SNP_5000.simseq.genome.fa
    ```

!!! terminal "code"

    ```
    #indexing the reference sample
    bwa index GCF_000191525.1_ASM19152v1_genomic.fna
    ```

!!! success "Output"

    ```
    [bwa_index] Pack FASTA... 0.01 sec
    [bwa_index] Construct BWT for the packed sequence...
    [bwa_index] 0.39 seconds elapse.
    [bwa_index] Update BWT... 0.02 sec
    [bwa_index] Pack forward-only FASTA... 0.00 sec
    [bwa_index] Construct SA from BWT and Occ... 0.24 sec
    [main] Version: 0.7.17-r1188
    [main] CMD: bwa index GCF_000191525.1_ASM19152v1_genomic.fna
    [main] Real time: 1.131 sec; CPU: 0.673 sec
    ```

In order to simulate a real sequencing experiment, we'll simulate the short reads too from the simulated full sequence using `wgsim` and map those reads to the reference sequence using `bwa`. Since the length of the each sequence is around 2.3 million, 0.7 millions of 100bp reads will give 30x read depth. (700000x100/2300000 ~ 30)

<!--
Why can't the following be a loop?

!!! terminal "code"

```bash
for simseq in *simseq.genome.fa; do
  # Get base name
  name=$(basename ${simseq} .simseq.genome.fa)
  # Simulate reads
  wgsim -N675000 -1100 -2100 ${simseq} ${name}.read1.fq ${name}.read2.fq
  # Specify read group header line
  rg_head="@RG\tID:${name}\tSM:${name}\tLB:L1"
  # Map simulated reads to reference genome
  bwa mem -R ${rg_head} \
    GCF_000191525.1_ASM19152v1_genomic.fna \
    ${name}.read1.fq ${name}.read2.fq \
    > ${name}.sam
  # Convert SAM to sorted BAM
  samtools view -bS ${name}.sam \
    | samtools sort - \
    > ${name}.bam
  # Call variants
  bcftools mpileup -Ou -f GCF_000191525.1_ASM19152v1_genomic.fna ${name}.bam \
    | bcftools call -vm -Oz -o ${name}.bwa.30x.100R.vcf.gz
  # Index variants
  bcftools index ${name}.bwa.30x.100R.vcf.gz
done
```
-->



!!! terminal "code"

    ```bash
    # Creating VCF files for each sample file considering GCF_000191525.1_ASM19152v1_genomic.fna as reference and simulating 30x 100bp reads. 
    #export OMP_NUM_THREADS=1
    wgsim -N675000 -1100 -2100 Simulation_SNP_5000.simseq.genome.fa Simulation_SNP_5000.read1.fq Simulation_SNP_5000.read2.fq 
    bwa mem -R "@RG\tID:Simulation_SNP_5000\tSM:Simulation_SNP_5000\tLB:L1" GCF_000191525.1_ASM19152v1_genomic.fna Simulation_SNP_5000.read1.fq Simulation_SNP_5000.read2.fq > Simulation_SNP_5000.sam
    samtools view -bS Simulation_SNP_5000.sam | samtools sort - > Simulation_SNP_5000.bam
    bcftools mpileup -Ou -f GCF_000191525.1_ASM19152v1_genomic.fna Simulation_SNP_5000.bam | bcftools call -vmO z -o Simulation_SNP_5000.bwa.30x.100R.vcf.gz
    bcftools index Simulation_SNP_5000.bwa.30x.100R.vcf.gz 
    ```

!!! terminal "code"

    ```bash    
    wgsim -N675000 -1100 -2100 Simulation_INDEL_5000.simseq.genome.fa Simulation_INDEL_5000.read1.fq Simulation_INDEL_5000.read2.fq 
    bwa mem -R "@RG\tID:Simulation_INDEL_5000\tSM:Simulation_INDEL_5000\tLB:L1" GCF_000191525.1_ASM19152v1_genomic.fna Simulation_INDEL_5000.read1.fq Simulation_INDEL_5000.read2.fq > Simulation_INDEL_5000.sam
    samtools view -bS Simulation_INDEL_5000.sam | samtools sort - > Simulation_INDEL_5000.bam
    bcftools mpileup -Ou -f GCF_000191525.1_ASM19152v1_genomic.fna Simulation_INDEL_5000.bam | bcftools call -vmO z -o Simulation_INDEL_5000.bwa.30x.100R.vcf.gz
    bcftools index Simulation_INDEL_5000.bwa.30x.100R.vcf.gz 
    ```

!!! terminal "code"
    ```bash 
    wgsim -N675000 -1100 -2100 Simulation_SNP_4000_INDEL_4000.simseq.genome.fa Simulation_SNP_4000_INDEL_4000.read1.fq Simulation_SNP_4000_INDEL_4000.read2.fq 
    bwa mem -R "@RG\tID:Simulation_SNP_4000_INDEL_4000\tSM:Simulation_SNP_4000_INDEL_4000\tLB:L1" GCF_000191525.1_ASM19152v1_genomic.fna Simulation_SNP_4000_INDEL_4000.read1.fq  Simulation_SNP_4000_INDEL_4000.read2.fq > Simulation_SNP_4000_INDEL_4000.sam
    samtools view -bS Simulation_SNP_4000_INDEL_4000.sam | samtools sort - > Simulation_SNP_4000_INDEL_4000.bam
    bcftools mpileup -Ou -f GCF_000191525.1_ASM19152v1_genomic.fna Simulation_SNP_4000_INDEL_4000.bam | bcftools call -vmO z -o Simulation_SNP_4000_INDEL_4000.bwa.30x.100R.vcf.gz
    bcftools index Simulation_SNP_4000_INDEL_4000.bwa.30x.100R.vcf.gz
    ```

!!! terminal "code"
    ```bash    
    wgsim -N675000 -1100 -2100 Simulation_SNP_4000_INDEL_4000_INV_4.simseq.genome.fa Simulation_SNP_4000_INDEL_4000_INV_4.read1.fq Simulation_SNP_4000_INDEL_4000_INV_4.read2.fq
    bwa mem -R "@RG\tID:Simulation_SNP_4000_INDEL_4000_INV_4\tSM:Simulation_SNP_4000_INDEL_4000_INV_4\tLB:L1" GCF_000191525.1_ASM19152v1_genomic.fna Simulation_SNP_4000_INDEL_4000_INV_4.read1.fq Simulation_SNP_4000_INDEL_4000_INV_4.read2.fq > Simulation_SNP_4000_INDEL_4000_INV_4.sam
    samtools view -bS Simulation_SNP_4000_INDEL_4000_INV_4.sam | samtools sort - > Simulation_SNP_4000_INDEL_4000_INV_4.bam
    bcftools mpileup -Ou -f GCF_000191525.1_ASM19152v1_genomic.fna Simulation_SNP_4000_INDEL_4000_INV_4.bam | bcftools call -vmO z -o Simulation_SNP_4000_INDEL_4000_INV_4.bwa.30x.100R.vcf.gz
    bcftools index Simulation_SNP_4000_INDEL_4000_INV_4.bwa.30x.100R.vcf.gz
    ```

!!! terminal "code"
    ```bash    
    wgsim -N675000 -1100 -2100 Simulation_SNP_4000_INDEL_4000_CNV_4.simseq.genome.fa Simulation_SNP_4000_INDEL_4000_CNV_4.read1.fq Simulation_SNP_4000_INDEL_4000_CNV_4.read2.fq
    bwa mem -R "@RG\tID:Simulation_SNP_4000_INDEL_4000_CNV_4\tSM:Simulation_SNP_4000_INDEL_4000_CNV_4\tLB:L1" GCF_000191525.1_ASM19152v1_genomic.fna Simulation_SNP_4000_INDEL_4000_CNV_4.read1.fq Simulation_SNP_4000_INDEL_4000_CNV_4.read2.fq > Simulation_SNP_4000_INDEL_4000_CNV_4.sam
    samtools view -bS Simulation_SNP_4000_INDEL_4000_CNV_4.sam | samtools sort - > Simulation_SNP_4000_INDEL_4000_CNV_4.bam
    bcftools mpileup -Ou -f GCF_000191525.1_ASM19152v1_genomic.fna Simulation_SNP_4000_INDEL_4000_CNV_4.bam | bcftools call -vmO z -o Simulation_SNP_4000_INDEL_4000_CNV_4.bwa.30x.100R.vcf.gz
    bcftools index Simulation_SNP_4000_INDEL_4000_CNV_4.bwa.30x.100R.vcf.gz
    ```

<!--Update to become an array job? Turn this into admonition? Or assume learners know how to run this in parallel? -->
(You can find above script as a Slurm job script here [vc_bwa_compare.sh](https://github.com/nuzla/Pangenome-Graphs-Workshop/blob/main/Scripts/vc_bwa_compare.sh))

Now we can find the variant call stats using `bcftools stats`. 

!!! terminal "code"

    ```bash
    bcftools stats Simulation_SNP_5000.bwa.30x.100R.vcf.gz | head -30
    ```

    ??? success "Output"
        ```bash
        # This file was produced by bcftools stats (1.9+htslib-1.9) and can be plotted using plot-vcfstats.
        # The command line was: bcftools stats  Simulation_SNP_5000.bwa.30x.100R.vcf.gz
        #
        # Definition of sets:
        # ID    [2]id   [3]tab-separated file names
        ID      0       Simulation_SNP_5000.bwa.30x.100R.vcf.gz
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
        SN      0       number of records:      6740
        SN      0       number of no-ALTs:      0
        SN      0       number of SNPs: 6395
        SN      0       number of MNPs: 0
        SN      0       number of indels:       345
        SN      0       number of others:       0
        SN      0       number of multiallelic sites:   1
        ```

<!-- Results can differ, either set seed for `wgsim` for consistency or explain stochastic nature of simulation! -->

# Mapping Reads using `vg giraffe` (Graph Method)
Here we map reads to a pangenome graph instead of single linear reference sequence. For example we'll consider the first case `Simulation_INDEL_5000.simseq.genome.fa`. We can build a graph considering the reference sequence `GCF_000191525.1_ASM19152v1_genomic.fna` and the ground truth VCF file 

!!! terminal "code"

    ```bash
    #Load additional model need for vg
    module load vg/1.46.0
    
    #Copy the vcf files into the folder
    cp ../Simulation_*.vcf ./
    ls -1trhs Simulation_*.vcf
    ```

!!! success "Output"

    ```
    768K Simulation_SNP_5000.refseq2simseq.SNP.vcf
    768K Simulation_INDEL_5000.refseq2simseq.INDEL.vcf
    768K Simulation_SNP_4000_INDEL_4000.refseq2simseq.SNP.vcf
    768K Simulation_SNP_4000_INDEL_4000.refseq2simseq.INDEL.vcf
     512 Simulation_SNP_4000_INDEL_4000_INV_4.refseq2simseq.inversion.vcf
    256K Simulation_SNP_4000_INDEL_4000_CNV_4.refseq2simseq.CNV.vcf
    ```

<!-- Loop/array job? -->
!!! terminal "code"

    ```bash
    #create tabix index
    bgzip Simulation_SNP_5000.refseq2simseq.SNP.vcf
    tabix Simulation_SNP_5000.refseq2simseq.SNP.vcf.gz
    
    #create the graph and index (-p is prefix for filenames)
    vg autoindex --workflow giraffe -r GCF_000191525.1_ASM19152v1_genomic.fna -v Simulation_SNP_5000.refseq2simseq.SNP.vcf.gz -p VG_SNP_5000
    
    #Map reads to the graph
    vg giraffe -Z VG_SNP_5000.giraffe.gbz -f Simulation_INDEL_5000.read1.fq -f Simulation_INDEL_5000.read2.fq -o SAM > Simulation_VG_SNP_5000.sam
    
    #Follow same procedure for VCF file creation
    samtools view -bS Simulation_VG_SNP_5000.sam | samtools sort - > Simulation_VG_SNP_5000.bam
    bcftools mpileup -Ou -f GCF_000191525.1_ASM19152v1_genomic.fna Simulation_VG_SNP_5000.bam | bcftools call -vmO z -o Simulation_VG_SNP_5000.giraffe.30x.100R.vcf.gz
    bcftools index Simulation_VG_SNP_5000.giraffe.30x.100R.vcf.gz 
    ```
    
We can follow the same procedure for the rest of the samples and generate VCF files. 

