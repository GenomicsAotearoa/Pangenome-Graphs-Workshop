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
    cp ./5NM_2K94/5NM*.gfa ./vg_deconstruct/5NM_2k94.gfa
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
    tail -5 5NM_2k94.gfa | less -S
    ```
    !!! success "Output"
        
        ```
        P       NC_003112.2     85316+,85318+,85319+,85321+,85322+,85323+,85325+,85327+,85328+,85330+,85331+,85333+,85334+,85336+,85337+,85
        P       NC_017518.1     85316+,85317+,85319+,85320+,85322+,85323+,85325+,85326+,85328+,85329+,85331+,85332+,85334+,85335+,85337+,85
        P       NZ_CP007668.1   1+,4+,5+,6+,8+,9+,11+,12+,14+,15+,17+,18+,20+,21+,23+,25+,26+,27+,29+,31+,32+,34+,35+,37+,38+,39+,41+,43+,4
        P       NZ_CP016880.1   2+,133478+,133479+,133481+,133482+,133483+,133485+,133486+,133488+,133489+,133490+,133492+,133493+,133495+,
        P       NZ_CP020423.2   3+,175915+,175916+,175918+,175919+,175921+,175922+,175924+,175925+,175926+,175928+,175929+,175931+,175932+, 81
        ```


!!! terminal "code"

    ```bash
    #use vg deconstruct the graph into VCF based on the first path NC_003112.2
    #-e, --path-traversals    Only consider traversals that correspond to paths in the graph.
    #-a, --all-snarls         Process all snarls, including nested snarls (by default only top-level snarls reported).
    
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

??? infor "The following is a SLURM script to deconstruct graph into vcf files"  

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
??? infor "submit the script using the `sbatch` command as follows. Take note of the job ID for tracking the run."

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
    #check the vcf 
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
### check the complex variation in vcf files
!!! terminal "code"

    ```bash
    #check complex variation  
    awk 'length($4) > 2' 5NM_2k94aep1.vcf |head -100 |less -S 
    ```
??? success "Output"
        
        ```
        ##INFO=<ID=AC,Number=A,Type=Integer,Description="Total number of alternate alleles in called genotypes">
        ##INFO=<ID=NS,Number=1,Type=Integer,Description="Number of samples with data">
        ##INFO=<ID=AN,Number=1,Type=Integer,Description="Total number of alleles in called genotypes">
        ##INFO=<ID=LV,Number=1,Type=Integer,Description="Level in the snarl tree (0=top level)">
        ##INFO=<ID=PS,Number=1,Type=String,Description="ID of variant corresponding to parent snarl">
        ##INFO=<ID=AT,Number=R,Type=String,Description="Allele Traversal as path in graph">
        #CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT  NC_017518.1     NZ_CP007668.1   NZ_CP016880.1   NZ_CP020423.2
        NC_003112.2     1166    >85479>85482    TCT     TG      60      .       AC=4;AF=1;AN=4;AT=>85479>85481>85482,>85479>85480>85482;NS=4;LV=0       GT      1       1       1       1
        NC_003112.2     1188    >85490>85493    CGT     GAC     60      .       AC=4;AF=1;AN=4;AT=>85490>85492>85493,>85490>85491>85493;NS=4;LV=0       GT      1       1       1       1
        NC_003112.2     1197    >85493>85496    GAC     TTT     60      .       AC=4;AF=1;AN=4;AT=>85493>85495>85496,>85493>85494>85496;NS=4;LV=0       GT      1       1       1       1
        NC_003112.2     1288    >85548>85550    TCG     T       60      .       AC=2;AF=0.5;AN=4;AT=>85548>85549>85550,>85548>85550;NS=4;LV=0   GT      0       0       1       1
        NC_003112.2     1461    >85601>85604    TGTGAAGAATTCA   CGTAAAGAACTCG   60      .       AC=4;AF=1;AN=4;AT=>85601>85603>85604,>85601>85602>85604;NS=4;LV=0       GT      1       1       1       1
        NC_003112.2     2421    >85870>85879    CCTTGTTTTCAATGCTTCGGCACGCGGAACAGTGTATCACGCGCCGCCGACCGAATTCCTTCGGGATTGCGTCCAAAAA CT,CCTTGTTTTCAATGCTTCGGCACGCGGAACAGTGTATCACGCGCCGCCGACCGGATTTCTTCGGGATTGCGTCCAAAAA      60      .       AC=2,2
        NC_003112.2     2522    >85895>85898    TGA     CCC     60      .       AC=2;AF=0.5;AN=4;AT=>85895>85897>85898,>85895>85896>85898;NS=4;LV=0     GT      0       0       1       1
        NC_003112.2     2531    >85901>85904    CCC     CG      60      .       AC=2;AF=0.5;AN=4;AT=>85901>85903>85904,>85901>85902>85904;NS=4;LV=0     GT      0       0       1       1
        NC_003112.2     3725    >86007>86010    AAC     GGT     60      .       AC=2;AF=0.5;AN=4;AT=>86007>86008>86010,>86007>86009>86010;NS=4;LV=0     GT      0       0       1       1
        NC_003112.2     4158    >86067>86069    CTTTT   C       60      .       AC=4;AF=1;AN=4;AT=>86067>86068>86069,>86067>86069;NS=4;LV=0     GT      1       1       1       1
        NC_003112.2     4174    >86077>86081    CCGG    CA,CGG  60      .       AC=2,2;AF=0.5,0.5;AN=4;AT=>86077>86079>86080>86081,>86077>86078>86081,>86077>86080>86081;NS=4;LV=0      GT      2       2       1       1
        NC_003112.2     4497    >86125>86128    GGG     AAC     60      .       AC=1;AF=0.25;AN=4;AT=>86125>86126>86128,>86125>86127>86128;NS=4;LV=0    GT      0       0       1       0
        NC_003112.2     5425    >86315>86318    GTCT    CAAC    60      .       AC=2;AF=0.5;AN=4;AT=>86315>86316>86318,>86315>86317>86318;NS=4;LV=0     GT      1       1       0       0
        NC_003112.2     5713    >86384>86387    ATTGGTCGATGCCATT        CTTGGTCGATGCCATC        60      .       AC=3;AF=0.75;AN=4;AT=>86384>86385>86387,>86384>86386>86387;NS=4;LV=0    GT      1       1       1       0
        NC_003112.2     6244    >86433>86436    TTT     CAA     60      .       AC=4;AF=1;AN=4;AT=>86433>86435>86436,>86433>86434>86436;NS=4;LV=0       GT      1       1       1       1
        NC_003112.2     6266    >86436>86438    TTT     T       60      .       AC=4;AF=1;AN=4;AT=>86436>86437>86438,>86436>86438;NS=4;LV=0     GT      1       1       1       1
        NC_003112.2     6654    >86555>86557    ATGTCTTTCACCACAGCTTCCGCATCGGCGGCAA      A       60      .       AC=4;AF=1;AN=4;AT=>86555>86556>86557,>86555>86557;NS=4;LV=0     GT      1       1       1       1
        NC_003112.2     8852    >86833>86835    GTATAGTGAATTAACAAAAATCGGGACAAGGCGGCGAAGCCGCAGACAGTACAGATAGTACAGAACCGATTCACTTGGTGCTTCAGCACCTTAGAGAATCGTTCTCTTTGAGCTAAGGCGAGGCAATACCGTACTGGTTTTTGTTAATCCACTA      G       60      .       AC=3;A
        NC_003112.2     9305    >86883>86885    CCC     C       60      .       AC=4;AF=1;AN=4;AT=>86883>86884>86885,>86883>86885;NS=4;LV=0     GT      1       1       1       1
        NC_003112.2     9316    >86888>86891    CCC     ATT     60      .       AC=4;AF=1;AN=4;AT=>86888>86890>86891,>86888>86889>86891;NS=4;LV=0       GT      1       1       1       1
        NC_003112.2     9341    >86904>86906    TCCGC   T       60      .       AC=4;AF=1;AN=4;AT=>86904>86905>86906,>86904>86906;NS=4;LV=0     GT      1       1       1       1
        NC_003112.2     11685   >87208>87211    CGCG    ATTT    60      .       AC=2;AF=0.5;AN=4;AT=>87208>87209>87211,>87208>87210>87211;NS=4;LV=0     GT      1       1       0       0
        NC_003112.2     13080   >87342>87372    AAAATGCCCCGCACCGTCTTTGCCGACCCAGTCGCAACACGGTTCGCCCTGCGACGTTTTGGCGGCAATCGCCTGAAAAATCGGCTTGACCGCATCCCAA            GAAATGTCCCGCACCATCCCTGCCGACCCAATCGCAACACGGCTCACCTTGAGGGGTTTTCGCAGCAATCGCCTGGAAAATCGGTT
        NC_003112.2     13473   >87464>87469    ATTA    GCCC,ATCC       60      .       AC=2,2;AF=0.5,0.5;AN=4;AT=>87464>87466>87467>87469,>87464>87465>87468>87469,>87464>87466>87468>87469;NS=4;LV=0  GT      2       2       1       1
        NC_003112.2     13849   >87531>87538    CAG     TGG,CGA 60      .       AC=2,1;AF=0.5,0.25;AN=4;AT=>87531>87533>87535>87537>87538,>87531>87532>87534>87537>87538,>87531>87533>87534>87536>87538;NS=4;LV=0       GT      2       0
        NC_003112.2     14058   >87606>87613    GGTC    GAGA,AGTA       60      .       AC=2,1;AF=0.5,0.25;AN=4;AT=>87606>87607>87610>87612>87613,>87606>87607>87608>87611>87613,>87606>87609>87610>87611>87613;NS=4;LV=0       GT      2
        NC_003112.2     14113   >87643>87846    AATTTGTGTATAAGTGGTGGAAAAAATGAGATTTGCGGGTAAATCTCACAATATTTCAGTCAGATAACTTTGGATTGCTTGTGTATAAGTAAACTTTCGGATGGGGATACGTAACGGAAACCTGTACCGCGTCATTCCCACGAACCTACATTCCGTCATTCCCACGAAAGTGGGAATGATGAAATTTTGA
        NC_003112.2     14291   >87652>87654    ATGAAATTTTGAGTTTTAGGAATTTATCGGGAGCAACAGAAACCGCTCCGCCGTCATTCCCGCGCAGGCGGGAATCTAGAACGTAAAATCTAAAGAAACCGTGTTGTAACGGCAGACCGATGCCGTCATTCCCGCGCAGGCGGGAATCTAGACCATTGGACAGCGGCAATATTCAAAGATTATCTGAAAG
        NC_003112.2     14598   >87660>87663    AGAACA  GTCGGT  60      .       AC=2;AF=1;AN=2;AT=>87660>87662>87663,>87660>87661>87663;NS=2;LV=1;PS=>87643>87846       GT      .       .       1       1
        NC_003112.2     14635   >87690>87692    GAGA    G       60      .       AC=2;AF=1;AN=2;AT=>87690>87691>87692,>87690>87692;NS=2;LV=1;PS=>87643>87846     GT      .       .       1       1
        NC_003112.2     14797   >87719>87722    AGAATA  GTCGGT  60      .       AC=2;AF=1;AN=2;AT=>87719>87720>87722,>87719>87721>87722;NS=2;LV=1;PS=>87643>87846       GT      .       .       1       1
        NC_003112.2     14828   >87743>87746    AAG     TTT     60      .       AC=2;AF=1;AN=2;AT=>87743>87745>87746,>87743>87744>87746;NS=2;LV=1;PS=>87643>87846       GT      .       .       1       1
        NC_003112.2     14834   >87749>87751    GGGA    G       60      .       AC=2;AF=1;AN=2;AT=>87749>87750>87751,>87749>87751;NS=2;LV=1;PS=>87643>87846     GT      .       .       1       1
        NC_003112.2     15023   >87770>87776    CCC     C,CTGCC,CCGCC   60      .               AC=1,1,1;AF=0.333333,0.333333,0.333333;AN=3;AT=>87770>87774>87776,>87770>87776,>87770>87771>87773>87774>87776,>87770>87772>87773>87774>87776;NS=3;LV=1
        NC_003112.2     15060   >87784>87787    GGTG    GAA     60      .       AC=3;AF=1;AN=3;AT=>87784>87786>87787,>87784>87785>87787;NS=3;LV=1;PS=>87643>87846       GT      1       .       1       1
        NC_003112.2     15064   >87787>87790    CGG     CA      60      .       AC=3;AF=1;AN=3;AT=>87787>87789>87790,>87787>87788>87790;NS=3;LV=1;PS=>87643>87846       GT      1       .       1       1
        NC_003112.2     15074   >87795>87797    TCGGA   T       60      .       AC=3;AF=1;AN=3;AT=>87795>87796>87797,>87795>87797;NS=3;LV=1;PS=>87643>87846     GT      1       .       1       1
        NC_003112.2     16456   >87941>88094    TTTGGAATTTCAATGCCTCAAGAATTTATCGGAAAAAACCAAAACCCTTCCGCCGTCATTCCCACGAAAGTGGGAATCTAGAAATGAAAAGCAGCAGGCATTTATCGGAAATGACCGAAACTGAACGGACTGGATTCCCGCTTTTGCGGGAATGACGGCGACAGGGTTGCTGTTATAGTGGATGAACAAA
        NC_003112.2     16475   >87961>87963    AAGAAT  A       60      .       AC=3;AF=1;AN=3;AT=>87961>87962>87963,>87961>87963;NS=3;LV=1;PS=>87941>88094     GT      .       1       1       1
        NC_003112.2     16485   >87966>87968    CGGA    C       60      .       AC=3;AF=1;AN=3;AT=>87966>87967>87968,>87966>87968;NS=3;LV=1;PS=>87941>88094     GT      .       1       1       1
        NC_003112.2     16499   >87977>87979    ACCC    A       60      .       AC=3;AF=1;AN=3;AT=>87977>87978>87979,>87977>87979;NS=3;LV=1;PS=>87941>88094     GT      .       1       1       1
        NC_003112.2     16507   >87979>87981    GCCG    G       60      .       AC=3;AF=1;AN=3;AT=>87979>87980>87981,>87979>87981;NS=3;LV=1;PS=>87941>88094     GT      .       1       1       1
        NC_003112.2     16537   >87987>87990    AAA     TTT     60      .       AC=3;AF=1;AN=3;AT=>87987>87989>87990,>87987>87988>87990;NS=3;LV=1;PS=>87941>88094       GT      .       1       1       1

        ```
