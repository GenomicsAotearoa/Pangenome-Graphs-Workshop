# 6. Variant calling
!!! info ""

- To detect both small and large variants among paths from the pangenome graph, we utilized the Variation Graph (VG) toolkit to deconstruct these variants into VCF files.
- To decompose the graph into a VCF file, we need to choose one path as a reference for comparison with the others (any path can serve as this reference)
- In fact, during the pangenome graph construction process, when the parameter -V 'NC_017518.1:#' is activated, the output file includes the VCF based on the NC_017518.1 reference.

## `vg deconstruct` graph to get the variations in vcf 
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

## check the vcf files
### check the statistics of vcf files
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

### check the vcf files

!!! terminal "code"

    ```bash
    #check the vcf statistics  
    head -100 5NM_2k94aep1.vcf |less -s 
    ```
    ??? success "Output"
        
        ```
        ##fileformat=VCFv4.2
        ##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
        ##INFO=<ID=AC,Number=A,Type=Integer,Description="Total number of alternate alleles in called genotypes">
        ##INFO=<ID=AF,Number=A,Type=Float,Description="Estimated allele frequency in the range (0,1]">
        ##INFO=<ID=NS,Number=1,Type=Integer,Description="Number of samples with data">
        ##INFO=<ID=AN,Number=1,Type=Integer,Description="Total number of alleles in called genotypes">
        ##INFO=<ID=LV,Number=1,Type=Integer,Description="Level in the snarl tree (0=top level)">
        ##INFO=<ID=PS,Number=1,Type=String,Description="ID of variant corresponding to parent snarl">
        ##INFO=<ID=AT,Number=R,Type=String,Description="Allele Traversal as path in graph">
        ##contig=<ID=NC_003112.2,length=2272360>
        #CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT  NC_017518.1     NZ_CP007668.1   NZ_CP016880.1   NZ_CP020423.2
        NC_003112.2     21      >85316>85319    C       T       60      .       AC=4;AF=1;AN=4;AT=>85316>85318>85319,>85316>85317>85319;NS=4;LV=0       GT      1       1       1       1
        NC_003112.2     77      >85319>85322    A       G       60      .       AC=4;AF=1;AN=4;AT=>85319>85321>85322,>85319>85320>85322;NS=4;LV=0       GT      1       1       1       1
        NC_003112.2     166     >85322>85325    G       A       60      .       AC=1;AF=0.25;AN=4;AT=>85322>85323>85325,>85322>85324>85325;NS=4;LV=0    GT      0       0       1       0
        NC_003112.2     226     >85325>85328    C       A       60      .       AC=4;AF=1;AN=4;AT=>85325>85327>85328,>85325>85326>85328;NS=4;LV=0       GT      1       1       1       1
        NC_003112.2     245     >85328>85331    C       T       60      .       AC=4;AF=1;AN=4;AT=>85328>85330>85331,>85328>85329>85331;NS=4;LV=0       GT      1       1       1       1
        NC_003112.2     279     >85331>85334    C       T       60      .       AC=4;AF=1;AN=4;AT=>85331>85333>85334,>85331>85332>85334;NS=4;LV=0       GT      1       1       1       1
        NC_003112.2     296     >85334>85337    A       G       60      .       AC=3;AF=0.75;AN=4;AT=>85334>85336>85337,>85334>85335>85337;NS=4;LV=0    GT      1       1       0       1
        NC_003112.2     322     >85337>85340    T       A       60      .       AC=3;AF=0.75;AN=4;AT=>85337>85339>85340,>85337>85338>85340;NS=4;LV=0    GT      1       1       0       1
        NC_003112.2     481     >85340>85343    C       T       60      .       AC=1;AF=0.25;AN=4;AT=>85340>85342>85343,>85340>85341>85343;NS=4;LV=0    GT      0       0       0       1
        NC_003112.2     517     >85343>85346    C       T       60      .       AC=4;AF=1;AN=4;AT=>85343>85345>85346,>85343>85344>85346;NS=4;LV=0       GT      1       1       1       1
        NC_003112.2     555     >85346>85349    G       A       60      .       AC=2;AF=0.5;AN=4;AT=>85346>85348>85349,>85346>85347>85349;NS=4;LV=0     GT      1       1       0       0
        NC_003112.2     582     >85349>85352    AA      GC      60      .       AC=4;AF=1;AN=4;AT=>85349>85351>85352,>85349>85350>85352;NS=4;LV=0       GT      1       1       1       1
        NC_003112.2     616     >85352>85355    T       C       60      .       AC=1;AF=0.25;AN=4;AT=>85352>85354>85355,>85352>85353>85355;NS=4;LV=0    GT      0       0       0       1
        NC_003112.2     660     >85355>85358    T       C       60      .       AC=1;AF=0.25;AN=4;AT=>85355>85356>85358,>85355>85357>85358;NS=4;LV=0    GT      0       0       1       0
        NC_003112.2     667     >85358>85361    A       C       60      .       AC=1;AF=0.25;AN=4;AT=>85358>85359>85361,>85358>85360>85361;NS=4;LV=0    GT      0       0       1       0
        NC_003112.2     670     >85361>85364    C       T       60      .       AC=1;AF=0.25;AN=4;AT=>85361>85363>85364,>85361>85362>85364;NS=4;LV=0    GT      0       0       1       0
        NC_003112.2     672     >85364>85367    G       T       60      .       AC=1;AF=0.25;AN=4;AT=>85364>85365>85367,>85364>85366>85367;NS=4;LV=0    GT      0       0       1       0
        NC_003112.2     730     >85367>85370    C       T       60      .       AC=4;AF=1;AN=4;AT=>85367>85369>85370,>85367>85368>85370;NS=4;LV=0       GT      1       1       1       1
        NC_003112.2     738     >85370>85373    A       G       60      .       AC=4;AF=1;AN=4;AT=>85370>85372>85373,>85370>85371>85373;NS=4;LV=0       GT      1       1       1       1
        NC_003112.2     745     >85373>85376    G       A       60      .       AC=4;AF=1;AN=4;AT=>85373>85375>85376,>85373>85374>85376;NS=4;LV=0       GT      1       1       1       1
        NC_003112.2     748     >85376>85379    G       A       60      .       AC=4;AF=1;AN=4;AT=>85376>85378>85379,>85376>85377>85379;NS=4;LV=0       GT      1       1       1       1
        NC_003112.2     751     >85379>85382    A       G       60      .       AC=3;AF=0.75;AN=4;AT=>85379>85381>85382,>85379>85380>85382;NS=4;LV=0    GT      1       1       0       1
        NC_003112.2     781     >85382>85385    A       C       60      .       AC=4;AF=1;AN=4;AT=>85382>85384>85385,>85382>85383>85385;NS=4;LV=0       GT      1       1       1       1
        NC_003112.2     790     >85385>85389    T       C,G     60      .       AC=3,1;AF=0.75,0.25;AN=4;AT=>85385>85388>85389,>85385>85386>85389,>85385>85387>85389;NS=4;LV=0  GT      1       1       2       1
        NC_003112.2     816     >85389>85392    C       A       60      .       AC=3;AF=0.75;AN=4;AT=>85389>85391>85392,>85389>85390>85392;NS=4;LV=0    GT      1       1       0       1
        NC_003112.2     829     >85392>85395    A       C       60      .       AC=3;AF=0.75;AN=4;AT=>85392>85394>85395,>85392>85393>85395;NS=4;LV=0    GT      1       1       0       1
        NC_003112.2     844     >85395>85398    T       G       60      .       AC=3;AF=0.75;AN=4;AT=>85395>85397>85398,>85395>85396>85398;NS=4;LV=0    GT      1       1       0       1
        NC_003112.2     847     >85398>85401    T       C       60      .       AC=3;AF=0.75;AN=4;AT=>85398>85400>85401,>85398>85399>85401;NS=4;LV=0    GT      1       1       0       1
        NC_003112.2     855     >85401>85404    A       G       60      .       AC=4;AF=1;AN=4;AT=>85401>85403>85404,>85401>85402>85404;NS=4;LV=0       GT      1       1       1       1
        NC_003112.2     859     >85404>85407    A       G       60      .       AC=3;AF=0.75;AN=4;AT=>85404>85406>85407,>85404>85405>85407;NS=4;LV=0    GT      1       1       0       1
        NC_003112.2     865     >85407>85410    T       C       60      .       AC=1;AF=0.25;AN=4;AT=>85407>85408>85410,>85407>85409>85410;NS=4;LV=0    GT      0       0       1       0
        NC_003112.2     877     >85410>85413    T       C       60      .       AC=3;AF=0.75;AN=4;AT=>85410>85412>85413,>85410>85411>85413;NS=4;LV=0    GT      1       1       0       1
        NC_003112.2     883     >85413>85416    G       A       60      .       AC=1;AF=0.25;AN=4;AT=>85413>85414>85416,>85413>85415>85416;NS=4;LV=0    GT      0       0       1       0
        NC_003112.2     899     >85416>85419    G       A       60      .       AC=1;AF=0.25;AN=4;AT=>85416>85417>85419,>85416>85418>85419;NS=4;LV=0    GT      0       0       1       0
        NC_003112.2     969     >85419>85422    A       G       60      .       AC=1;AF=0.25;AN=4;AT=>85419>85420>85422,>85419>85421>85422;NS=4;LV=0    GT      0       0       1       0
        NC_003112.2     981     >85422>85425    G       T       60      .       AC=3;AF=0.75;AN=4;AT=>85422>85424>85425,>85422>85423>85425;NS=4;LV=0    GT      1       1       1       0
        NC_003112.2     987     >85425>85428    T       C       60      .       AC=1;AF=0.25;AN=4;AT=>85425>85427>85428,>85425>85426>85428;NS=4;LV=0    GT      0       0       1       0
        NC_003112.2     996     >85428>85431    A       G       60      .       AC=1;AF=0.25;AN=4;AT=>85428>85430>85431,>85428>85429>85431;NS=4;LV=0    GT      0       0       1       0
        NC_003112.2     1005    >85431>85434    G       A       60      .       AC=1;AF=0.25;AN=4;AT=>85431>85432>85434,>85431>85433>85434;NS=4;LV=0    GT      0       0       1       0
        NC_003112.2     1008    >85434>85437    A       G       60      .       AC=3;AF=0.75;AN=4;AT=>85434>85436>85437,>85434>85435>85437;NS=4;LV=0    GT      1       1       1       0
        NC_003112.2     1016    >85437>85440    GG      AA      60      .       AC=1;AF=0.25;AN=4;AT=>85437>85438>85440,>85437>85439>85440;NS=4;LV=0    GT      0       0       1       0
        NC_003112.2     1023    >85440>85443    G       A       60      .       AC=2;AF=0.5;AN=4;AT=>85440>85442>85443,>85440>85441>85443;NS=4;LV=0     GT      1       1       0       0
        NC_003112.2     1032    >85443>85446    G       A       60      .       AC=2;AF=0.5;AN=4;AT=>85443>85445>85446,>85443>85444>85446;NS=4;LV=0     GT      1       1       0       0
        NC_003112.2     1041    >85446>85449    A       G       60      .       AC=2;AF=0.5;AN=4;AT=>85446>85447>85449,>85446>85448>85449;NS=4;LV=0     GT      1       1       0       0
        NC_003112.2     1044    >85449>85452    G       A       60      .       AC=1;AF=0.25;AN=4;AT=>85449>85451>85452,>85449>85450>85452;NS=4;LV=0    GT      0       0       1       0
        NC_003112.2     1056    >85452>85455    G       A       60      .       AC=1;AF=0.25;AN=4;AT=>85452>85454>85455,>85452>85453>85455;NS=4;LV=0    GT      0       0       1       0

        ```


