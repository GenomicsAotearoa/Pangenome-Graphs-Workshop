# 6. Variant calling
!!! info ""

- To detect both small and large variants among paths from the pangenome graph, we utilized the Variation Graph (VG) toolkit to deconstruct these variants into VCF files.
- To decompose the graph into a VCF file, we need to choose one path as a reference for comparison with the others (any path can serve as this reference)
- In fact, during the pangenome graph construction process, when the parameter -V 'NC_017518.1:#' is activated, the output file includes the VCF based on the NC_017518.1 reference.

### `vg deconstruct` graph to get the variations in vcf 
!!! info ""

Set up directory for VG and GFA files.

!!! terminal "code"

    ```bash
    # Return to working directory
    cd ~/pg_workshop

    # Create VG directory
    mkdir -p vg_deconstruct

    #copy the gfa file to the directory
    cp ./5NM_2k94/5NM*.gfa ./vg_deconstruct/5NM_2k94.fa
    cd ./vg_deconstruct
    ```
    
Load the necessary modules for an example run.
!!! terminal "code"

    ```bash
    module purge
    module load vg/1.46.0
    module load BCFtools/1.15.1-GCC-11.3.0
    ```

An example run to obtain VCF files from GFA.

!!! terminal "code"

    ```bash
    #check the paths in the graph using tail, which depends on the number of genomes. We have five input genomes for the 5NM dataset. 
    tail -5 5NM_2k94.gfa
    ```
    !!! success "Output"
        
        ```
        P       NC_003112.2     85316+,85318+,85319+,85321+,85322+,85323+,85325+,85327+,85328+,85330+,85331+,85333+,85334+,85336+,85337+,85
        P       NC_017518.1     85316+,85317+,85319+,85320+,85322+,85323+,85325+,85326+,85328+,85329+,85331+,85332+,85334+,85335+,85337+,85
        P       NZ_CP007668.1   1+,4+,5+,6+,8+,9+,11+,12+,14+,15+,17+,18+,20+,21+,23+,25+,26+,27+,29+,31+,32+,34+,35+,37+,38+,39+,41+,43+,4
        P       NZ_CP016880.1   2+,133478+,133479+,133481+,133482+,133483+,133485+,133486+,133488+,133489+,133490+,133492+,133493+,133495+,
        P       NZ_CP020423.2   3+,175915+,175916+,175918+,175919+,175921+,175922+,175924+,175925+,175926+,175928+,175929+,175931+,175932+, 81
        ```


??? terminal "code"

    ```bash
    #use vg deconstruct the graph into VCF based on the first path NC_003112.2
    vg deconstruct -p NC_003112.2 -a -e ./5NM_2k94.gfa > 5NM_2k94aep1.vcf
    #use bcftools stats to check the statistics for the vcf file 
    bcftools stats 5NM_2k94aep1.vcf > 5NM_2k94aep1.vcf_stats
    ```

    ```bash
    #use vg deconstruct the graph into VCF based on the second path NC_017518.1
    vg deconstruct -p NC_017518.1 -a -e ./5NM_2k94.gfa > 5NM_2k94aep2.vcf
    
    #use bcftools stats to check the statistics for the vcf file 
    bcftools stats 5NM_2k94aep2.vcf > 5NM_2k94aep2.vcf_stats
    ```

!!! terminal "code"
The following is a SLURM script to deconstruct graph into vcf files  

    ```bash
    #!/usr/bin/bash

    #SBATCH --account       ga03793
    #SBATCH --job-name      deconstruct_gfa
    #SBATCH --cpus-per-task 8
    #SBATCH --mem           4G
    #SBATCH --time          1:00:00

    module purge
    module load vg/1.46.0
    module load BCFtools/1.15.1-GCC-11.3.0


    #use vg deconstruct the graph into VCF based on the first path NC_003112.2
    vg deconstruct -p NC_003112.2 -a -e ./5NM_2k94.gfa > 5NM_2k94aep1.vcf
    bcftools stats 5NM_2k94aep1.vcf > 5NM_2k94aep1.vcf_stats


    #use vg deconstruct the graph into VCF based on the second path NC_017518.1
    vg deconstruct -p NC_017518.1 -a -e ./5NM_2k94.gfa > 5NM_2k94aep2.vcf
    bcftools stats 5NM_2k94aep2.vcf > 5NM_2k94aep2.vcf_stats

    
    ```
submit the script using the `sbatch` command as follows. Take note of the job ID for tracking the run.

    ```bash
    sbatch vg_decon.sl
    ```

### check the vcf files

!!! terminal "code"

    ```bash
    #check the vcf statistics  
    less -S 5NM_2k94aep1.vcf_stats
    ```
    ??? success "Output"
        
        ```
        # This file was produced by bcftools stats (1.15.1+htslib-1.15.1) and can be plotted using plot-vcfstats.
        # The command line was: bcftools stats  5NM_2k94aep1.vcf
        #
        # Definition of sets:
        # ID    [2]id   [3]tab-separated file names
        ID      0       5NM_2k94aep1.vcf
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
        SN      0       number of samples:      4
        SN      0       number of records:      75538
        SN      0       number of no-ALTs:      0
        SN      0       number of SNPs: 66088
        SN      0       number of MNPs: 6948
        SN      0       number of indels:       3284
        SN      0       number of others:       968
        SN      0       number of multiallelic sites:   3705
        SN      0       number of multiallelic SNP sites:       1344
        # TSTV, transitions/transversions:
        # TSTV  [2]id   [3]ts   [4]tv   [5]ts/tv        [6]ts (1st ALT) [7]tv (1st ALT) [8]ts/tv (1st ALT)
        ```


!!! terminal "code"

    ```bash
    #check the vcf statistics  
    less -S 5NM_2k94aep2.vcf_stats
    ```
    ??? success "Output"
        
        ```
        #only show the number of each type of variations
        SN      0       number of samples:      4
        SN      0       number of records:      76867
        SN      0       number of no-ALTs:      0
        SN      0       number of SNPs: 67338
        SN      0       number of MNPs: 7119
        SN      0       number of indels:       3311
        SN      0       number of others:       957
        SN      0       number of multiallelic sites:   3758
        SN      0       number of multiallelic SNP sites:       1335

        ```




