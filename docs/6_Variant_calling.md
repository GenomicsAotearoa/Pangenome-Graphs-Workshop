# 6. Variant calling
!!! info ""

- To detect both small and large variants among paths from the pangenome graph, we utilized the Variation Graph (VG) toolkit to deconstruct these variants into VCF files.
- To decompose the graph into a VCF file, we need to choose one path as a reference for comparison with the others (any path can serve as this reference)
- In fact, during the pangenome graph construction process, when the parameter -V 'NC_017518.1:#' is activated, the output file includes the VCF based on the NC_017518.1 reference.

## `vg deconstruct` graph to get the variations in vcf 
!!! info ""

Set up directory for VG and GFA files.

!!! terminal "code"
    
    - Return to working directory
    ```bash
    cd ~/pg_workshop
    ```

    - Create VG directory
    ```bash
    mkdir -p vg_deconstruct
    ```

    - Copy the gfa file to the directory
    ```bash
    cp ./5NM_2Kb94/5NM*.gfa ./vg_deconstruct/5NM_2Kb94.gfa
    cp ./5NM_2Kb94_k35/5NM*.gfa ././vg_deconstruct/5NM_2Kb94_k35.gfa
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

    - check the paths in the graph using tail, which depends on the number of genomes. We have five input genomes for the 5NM dataset. 
    ```bash
    tail -5 5NM_2Kb94.gfa | less -Sail -5 5NM_2k94.gfa | less -S
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
    #-H, --path-sep SEP       Obtain alt paths from the set of paths, assuming a path name hierarchy (e.g. SEP='#' and sample#phase#contig)
    ```
    
    - vg deconstruct for the 5NM_2Kb94.gfa using the path NC_003112.2 as reference 
    ```bash
    vg deconstruct -p NC_003112.2 -a -e -H AAAA ./5NM_2Kb94.gfa > 5NM_2Kb94aep1.vcf
    ```

    - use bcftools stats to check the statistics for the vcf file 
    ```
    bcftools stats 5NM_2Kb94aep1.vcf > 5NM_2Kb94aep1.vcf_stats
    ```
    - vg deconstruct for the 5NM_2Kb94_k35.gfa
    ```bash
    vg deconstruct -p NC_003112.2 -a -e -H AAAA ./5NM_2Kb94_k35.gfa > 5NM_2Kb94_k35aep1.vcf
    ```

    - use bcftools stats to check the statistics for the vcf file 
    ```
    bcftools stats 5NM_2Kb94_k35aep1.vcf > 5NM_2Kb94_k35aep1.vcf_stats    
    ```

    - use vg deconstruct the graph into VCF based on the second path NC_017518.1
    ```
    vg deconstruct -p NC_017518.1 -a -e -H AAAA ./5NM_2Kb94.gfa > 5NM_2Kb94aep2.vcf
    ```
    
    - use bcftools stats to check the statistics for the vcf file 
    ```bash
    bcftools stats 5NM_2Kb94aep2.vcf > 5NM_2Kb94aep2.vcf_stats
    ```


## check the vcf files

!!! terminal-2 "check the statistics of vcf files"

    ```bash
    less -S 5NM_2Kb94aep1.vcf_stats
    ```
    ??? success "Output"
        
        ```
        # This file was produced by bcftools stats (1.15.1+htslib-1.15.1) and can be plotted using plot-vcfstats.
        # The command line was: bcftools stats  5NM_2Kb94aep1.vcf
        #
        # Definition of s  ets:
        # ID    [2]id   [3]tab-separated file names
        ID      0       5NM_2Kb94aep1.vcf
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
        SN      0       number of records:      76185
        SN      0       number of no-ALTs:      0
        SN      0       number of SNPs: 66461
        SN      0       number of MNPs: 7044
        SN      0       number of indels:       3458
        SN      0       number of others:       995
        SN      0       number of multiallelic sites:   3735
        SN      0       number of multiallelic SNP sites:       1350
        ```
<br>

!!! terminal-2 "check the statistics of vcf files"

    ```bash
    less -S 5NM_2Kb94_k35aep1.vcf_stats
    ```
    ??? success "Output"
        
        ```
        #only show the number of each type of variations
        SN      0       number of samples:      4
        SN      0       number of records:      74293
        SN      0       number of no-ALTs:      0
        SN      0       number of SNPs: 64666
        SN      0       number of MNPs: 6963
        SN      0       number of indels:       3428
        SN      0       number of others:       1031
        SN      0       number of multiallelic sites:   3798
        SN      0       number of multiallelic SNP sites:       1290

        ```
<br>

!!! terminal-2 "Inspect the vcf files"

    ```bash
    head -100 5NM_2Kb94aep1.vcf |less -s 
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
        NC_003112.2     152     >1>4    C       T       60      .       AC=2;AF=0.5;AN=4;AT=>1>3>4,>1>2>4;NS=4;LV=0     GT      1       1       0       0
        NC_003112.2     510     >4>7    A       G       60      .       AC=2;AF=0.5;AN=4;AT=>4>6>7,>4>5>7;NS=4;LV=0     GT      0       0       1       1
        NC_003112.2     558     >7>10   A       G       60      .       AC=2;AF=0.5;AN=4;AT=>7>9>10,>7>8>10;NS=4;LV=0   GT      0       0       1       1
        NC_003112.2     954     >10>13  G       A       60      .       AC=2;AF=0.5;AN=4;AT=>10>11>13,>10>12>13;NS=4;LV=0       GT      1       1       0       0
        NC_003112.2     1139    >13>16  A       G       60      .       AC=1;AF=0.25;AN=4;AT=>13>14>16,>13>15>16;NS=4;LV=0      GT      0       0       1       0
        NC_003112.2     1411    >16>19  G       A       60      .       AC=2;AF=0.5;AN=4;AT=>16>18>19,>16>17>19;NS=4;LV=0       GT      1       1       0       0
        NC_003112.2     1539    >19>22  T       C       60      .       AC=2;AF=0.5;AN=4;AT=>19>21>22,>19>20>22;NS=4;LV=0       GT      1       1       0       0
        NC_003112.2     1561    >22>25  A       G       60      .       AC=4;AF=1;AN=4;AT=>22>23>25,>22>24>25;NS=4;LV=0 GT      1       1       1       1
        NC_003112.2     1630    >25>28  G       C       60      .       AC=4;AF=1;AN=4;AT=>25>26>28,>25>27>28;NS=4;LV=0 GT      1       1       1       1
        NC_003112.2     1674    >28>31  T       C       60      .       AC=4;AF=1;AN=4;AT=>28>29>31,>28>30>31;NS=4;LV=0 GT      1       1       1       1
        NC_003112.2     1781    >31>34  A       G       60      .       AC=4;AF=1;AN=4;AT=>31>32>34,>31>33>34;NS=4;LV=0 GT      1       1       1       1
        NC_003112.2     1809    >34>36  AT      A       60      .       AC=4;AF=1;AN=4;AT=>34>35>36,>34>36;NS=4;LV=0    GT      1       1       1       1
        NC_003112.2     1819    >36>38  AT      A       60      .       AC=4;AF=1;AN=4;AT=>36>37>38,>36>38;NS=4;LV=0    GT      1       1       1       1
        NC_003112.2     1894    >38>41  G       C       60      .       AC=4;AF=1;AN=4;AT=>38>39>41,>38>40>41;NS=4;LV=0 GT      1       1       1       1
        NC_003112.2     1909    >41>44  A       G       60      .       AC=4;AF=1;AN=4;AT=>41>42>44,>41>43>44;NS=4;LV=0 GT      1       1       1       1
        NC_003112.2     1974    >44>47  A       G       60      .       AC=4;AF=1;AN=4;AT=>44>45>47,>44>46>47;NS=4;LV=0 GT      1       1       1       1
        NC_003112.2     2069    >47>50  A       G       60      .       AC=4;AF=1;AN=4;AT=>47>48>50,>47>49>50;NS=4;LV=0 GT      1       1       1       1
        NC_003112.2     2073    >50>53  T       C       60      .       AC=2;AF=0.5;AN=4;AT=>50>51>53,>50>52>53;NS=4;LV=0       GT      0       0       1       1
        NC_003112.2     2079    >53>56  T       C       60      .       AC=4;AF=1;AN=4;AT=>53>54>56,>53>55>56;NS=4;LV=0 GT      1       1       1       1
        NC_003112.2     2084    >56>59  AT      GG      60      .       AC=4;AF=1;AN=4;AT=>56>57>59,>56>58>59;NS=4;LV=0 GT      1       1       1       1
        NC_003112.2     2089    >59>62  T       C       60      .       AC=4;AF=1;AN=4;AT=>59>60>62,>59>61>62;NS=4;LV=0 GT      1       1       1       1
        NC_003112.2     2099    >62>64  A       AC      60      .       AC=4;AF=1;AN=4;AT=>62>64,>62>63>64;NS=4;LV=0    GT      1       1       1       1
        NC_003112.2     2178    >64>67  G       A       60      .       AC=4;AF=1;AN=4;AT=>64>65>67,>64>66>67;NS=4;LV=0 GT      1       1       1       1
        NC_003112.2     2200    >67>70  A       G       60      .       AC=4;AF=1;AN=4;AT=>67>68>70,>67>69>70;NS=4;LV=0 GT      1       1       1       1
        NC_003112.2     2290    >70>73  C       T       60      .       AC=4;AF=1;AN=4;AT=>70>71>73,>70>72>73;NS=4;LV=0 GT      1       1       1       1
        NC_003112.2     2359    >73>76  T       C       60      .       AC=4;AF=1;AN=4;AT=>73>74>76,>73>75>76;NS=4;LV=0 GT      1       1       1       1
        NC_003112.2     2362    >76>79  C       T       60      .       AC=4;AF=1;AN=4;AT=>76>77>79,>76>78>79;NS=4;LV=0 GT      1       1       1       1
        NC_003112.2     2392    >79>84  CG      AG,AA   60      .       AC=2,2;AF=0.5,0.5;AN=4;AT=>79>80>82>84,>79>81>82>84,>79>81>83>84;NS=4;LV=0      GT      2       2       1       1
        NC_003112.2     2491    >84>87  G       A       60      .       AC=2;AF=0.5;AN=4;AT=>84>85>87,>84>86>87;NS=4;LV=0       GT      0       0       1       1
        NC_003112.2     2494    >87>90  T       C       60      .       AC=2;AF=0.5;AN=4;AT=>87>88>90,>87>89>90;NS=4;LV=0       GT      0       0       1       1
        NC_003112.2     2503    >90>93  A       C       60      .       AC=4;AF=1;AN=4;AT=>90>91>93,>90>92>93;NS=4;LV=0 GT      1       1       1       1
        NC_003112.2     2507    >93>96  G       A       60      .       AC=4;AF=1;AN=4;AT=>93>94>96,>93>95>96;NS=4;LV=0 GT      1       1       1       1
        NC_003112.2     2515    >96>99  T       C       60      .       AC=4;AF=1;AN=4;AT=>96>97>99,>96>98>99;NS=4;LV=0 GT      1       1       1       1
        NC_003112.2     2520    >99>101 CA      C       60      .       AC=4;AF=1;AN=4;AT=>99>100>101,>99>101;NS=4;LV=0 GT      1       1       1       1
        NC_003112.2     2523    >101>106        CC      CAA,CAG 60      .       AC=3,1;AF=0.75,0.25;AN=4;AT=>101>102>106,>101>103>105>106,>101>103>104>106;NS=4;LV=0    GT      2       1       1       1

        ```
<br>

!!! terminal-2 "check the complex variation in vcf files"

    ```bash
    awk 'length($4) > 2' 5NM_2Kb94aep1.vcf |head -100 |less -S 
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
            NC_003112.2     3441    >250>600        CCAACCTCGCCAAAGTCCGCAAACAAGTAACTGCTTTGTGCAATAAATACCCCGTTTACGGCGCGTAAGCCTTTTTAAAAATATTCCGCCAAGCAATCCAATGCCGCCTGAAATCTCATAATGTTTCAGGCGGAAACCTTTGCAAAAATCCCCAAAATCCCCTAAATTCCCACCAAG
            NC_003112.2     4164    >392>397        GCG     AAG,AAA 60      .       AC=1,1;AF=0.5,0.5;AN=2;AT=>392>393>395>397,>392>394>395>397,>392>394>396>397;NS=2;LV=1;PS=>250>600      GT      2       1       .       .
            NC_003112.2     4178    >397>402        ATCAAGAAAAACGGC AT,ACAAAGAAAAACGGC      60      .       AC=1,1;AF=0.5,0.5;AN=2;AT=>397>398>399>401>402,>397>398>402,>397>400>401>402;NS=2;LV=1;PS=>250>600      GT      2       1
            NC_003112.2     4394    >426>431        TTG     TTA,CAG 60      .       AC=1,1;AF=0.5,0.5;AN=2;AT=>426>428>429>431,>426>428>430>431,>426>427>429>431;NS=2;LV=1;PS=>250>600      GT      2       1       .       .
            NC_003112.2     4918    >518>520        TCC     T       60      .       AC=1;AF=0.5;AN=2;AT=>518>519>520,>518>520;NS=2;LV=1;PS=>250>600 GT      1       0       .       .
            NC_003112.2     5116    >529>553        TCCCGTCATTCCCGCGCAGGCGGGAATCTAGGTTTGTCGGCACGGAAACTTATCGGGTAAAACGGTTTCTTTAGATTTTACGTTCTAGATTCCCGCCTGCGCGGGAATGACGATGAAAAGATTGTTGTCGCTTCGGATAAATTTTTGTCGCGTTGGGTTCTAGATTCCCGCCTGCGC
            NC_003112.2     5441    >553>556        TCCCC   T,TCCC  60      .       AC=1,1;AF=0.5,0.5;AN=2;AT=>553>554>555>556,>553>556,>553>555>556;NS=2;LV=1;PS=>250>600  GT      2       1       .       .
            NC_003112.2     6462    >600>603        ACT     AA      60      .       AC=3;AF=0.75;AN=4;AT=>600>602>603,>600>601>603;NS=4;LV=0        GT      1       0       1       1
            NC_003112.2     8056    >938>1174       AGCTGCGCCGCCAGCGTGCAGAACGCCGCAACGTACAGGGACAAGTCAACTTCAAGCTCGATAGCGGTGAGAAAAGTGGCAAAATCATCGCCGAATTGGAACACGCTTCGTTTGCCTATGGCGGCAAAGTCATTATGGACAAATTCTCCGCTATCTTGCAGCGCGGCGACAAAATCG
            NC_003112.2     40957   >1264>1267      CACCCAGTTGCGCCAAAGCTGCGCCATCCCGCTC      TACACAGTTACGCCAAAGTTGCGCCATCCCGCTT      60      .       AC=2;AF=0.5;AN=4;AT=>1264>1266>1267,>1264>1265>1267;NS=4;LV=0   GT      0       0
            NC_003112.2     42763   >1445>1450      TGG     CAA,CAG 60      .       AC=2,1;AF=0.5,0.25;AN=4;AT=>1445>1447>1449>1450,>1445>1446>1448>1450,>1445>1446>1449>1450;NS=4;LV=0     GT      2       0       1       1
            NC_003112.2     42797   >1465>1470      CGCG    CATA,TGCG       60      .       AC=2,1;AF=0.5,0.25;AN=4;AT=>1465>1466>1469>1470,>1465>1466>1467>1470,>1465>1468>1469>1470;NS=4;LV=0     GT      0       2       1       1
            NC_003112.2     43307   >1650>1653      GAG     TTC     60      .       AC=1;AF=0.25;AN=4;AT=>1650>1652>1653,>1650>1651>1653;NS=4;LV=0  GT      0       1       0       0
            NC_003112.2     45378   >1733>1738      GGCCCGCTTCAGGGGC        GGCGCGCTTCAGGGGC,G      60      .       AC=2,2;AF=0.5,0.5;AN=4;AT=>1733>1734>1736>1737>1738,>1733>1734>1735>1737>1738,>1733>1738;NS=4;LV=0      GT      2
            NC_003112.2     47380   >1813>1815      TCG     T       60      .       AC=4;AF=1;AN=4;AT=>1813>1814>1815,>1813>1815;NS=4;LV=0  GT      1       1       1       1
            NC_003112.2     47406   >1837>1839      AACGGAAT        A       60      .       AC=4;AF=1;AN=4;AT=>1837>1838>1839,>1837>1839;NS=4;LV=0  GT      1       1       1       1
            NC_003112.2     47431   >1851>1854      AAA     CTT     60      .       AC=4;AF=1;AN=4;AT=>1851>1852>1854,>1851>1853>1854;NS=4;LV=0     GT      1       1       1       1
            NC_003112.2     47447   >1861>1863      ACG     A       60      .       AC=4;AF=1;AN=4;AT=>1861>1862>1863,>1861>1863;NS=4;LV=0  GT      1       1       1       1
            NC_003112.2     47451   >1863>1866      TCCG    TT      60      .       AC=4;AF=1;AN=4;AT=>1863>1865>1866,>1863>1864>1866;NS=4;LV=0     GT      1       1       1       1
            NC_003112.2     47457   >1866>1869      GAA     TCC     60      .       AC=4;AF=1;AN=4;AT=>1866>1868>1869,>1866>1867>1869;NS=4;LV=0     GT      1       1       1       1
            NC_003112.2     47464   >1872>1874      TCCATTGCCGCATTGTCCAAGCAGTTTCCCTTGCGGGACATACTCTGAACCAGACCGTTGCCTTTCAACTGCTTTTG   T       60      .       AC=4;AF=1;AN=4;AT=>1872>1873>1874,>1872>1874;NS=4;LV=0  GT      1       1
            NC_003112.2     48557   >1944>1949      TCG     CGG,CGA 60      .       AC=2,1;AF=0.5,0.25;AN=4;AT=>1944>1946>1947>1949,>1944>1945>1947>1949,>1944>1945>1948>1949;NS=4;LV=0     GT      0       2       1       1
            NC_003112.2     48645   >1982>1985      CAC     AGA     60      .       AC=1;AF=0.25;AN=4;AT=>1982>1984>1985,>1982>1983>1985;NS=4;LV=0  GT      0       0       0       1
            NC_003112.2     48936   >2036>2040      TAT     TGT,TG  60      .       AC=2,1;AF=0.5,0.25;AN=4;AT=>2036>2037>2039>2040,>2036>2038>2039>2040,>2036>2038>2040;NS=4;LV=0  GT      0       1       2       1
            NC_003112.2     48942   >2040>2042      CCT     C       60      .       AC=1;AF=0.25;AN=4;AT=>2040>2041>2042,>2040>2042;NS=4;LV=0       GT      0       0       1       0
            NC_003112.2     49330   >2124>2129      TGT     TGG,ATT 60      .       AC=3,1;AF=0.75,0.25;AN=4;AT=>2124>2126>2127>2129,>2124>2126>2128>2129,>2124>2125>2127>2129;NS=4;LV=0    GT      1       2       1       1
            NC_003112.2     49436   >2170>2268      ACTTGTCTGACATGGAAAAATCCCTGTATTGAATTAAAAATCAATACAGGGATTGTAGGAAAGGCCGTCTGACTAAGCCTTTAATACGGGTTAAAACTTAATCAGTAGAGAGAGATGTGAGGATGATTTTTTTAGGCTTACGAGAGCCATTTTGCTTTAAGTCAAACTCAACTGTTA
            NC_003112.2     49437   >2173>2175      CTT     C       60      .       AC=3;AF=1;AN=3;AT=>2173>2174>2175,>2173>2175;NS=3;LV=1;PS=>2170>2268    GT      .       1       1       1
            NC_003112.2     50804   >2203>2205      AGAT    A       60      .       AC=1;AF=1;AN=1;AT=>2203>2204>2205,>2203>2205;NS=1;LV=1;PS=>2170>2268    GT      1       .       .       .
            NC_003112.2     50904   >2220>2223      ATTC    GGCA    60      .       AC=1;AF=1;AN=1;AT=>2220>2222>2223,>2220>2221>2223;NS=1;LV=1;PS=>2170>2268       GT      1       .       .       .
            NC_003112.2     53387   >2342>2345      CAA     ATG     60      .       AC=1;AF=0.25;AN=4;AT=>2342>2344>2345,>2342>2343>2345;NS=4;LV=0  GT      1       0       0       0
            NC_003112.2     53393   >2348>2351      ACAA    TTTC    60      .       AC=1;AF=0.25;AN=4;AT=>2348>2350>2351,>2348>2349>2351;NS=4;LV=0  GT      1       0       0       0
            NC_003112.2     53414   >2368>2370      CGAT    C       60      .       AC=1;AF=0.25;AN=4;AT=>2368>2369>2370,>2368>2370;NS=4;LV=0       GT      1       0       0       0
            NC_003112.2     53471   >2391>2394      AGCG    CCAC    60      .       AC=1;AF=0.25;AN=4;AT=>2391>2392>2394,>2391>2393>2394;NS=4;LV=0  GT      1       0       0       0
            NC_003112.2     53493   >2403>2406      TCAA    CTTG    60      .       AC=1;AF=0.25;AN=4;AT=>2403>2404>2406,>2403>2405>2406;NS=4;LV=0  GT      1       0       0       0
            NC_003112.2     53604   >2454>2456      TGGGCTG T       60      .       AC=1;AF=0.25;AN=4;AT=>2454>2455>2456,>2454>2456;NS=4;LV=0       GT      1       0       0       0
            NC_003112.2     53746   >2505>2508      CTGGTGGGTTTGGCCGAAGGG   TTGGTGGGACTTGCCGAAGGA   60      .       AC=3;AF=0.75;AN=4;AT=>2505>2507>2508,>2505>2506>2508;NS=4;LV=0  GT      1       1       1       0
            NC_003112.2     55432   >2590>2597      CACG    CTC,CCCGAACAACCCG       60      .       AC=2,2;AF=0.5,0.5;AN=4;AT=>2590>2592>2593>2596>2597,>2590>2591>2594>2597,>2590>2593>2594>2595>2596>2597;NS=4;LV=0       GT      2
            NC_003112.2     55984   >2637>2639      AAAACCACAACC    A       60      .       AC=3;AF=0.75;AN=4;AT=>2637>2638>2639,>2637>2639;NS=4;LV=0       GT      1       0       1       1
    
            ```

??? terminal-2 "bcftools `isec` to check the difference for 5NM based on two settings, -S 2000 -p 94 -k 19 and -S 2000 -p 94 -k 35"
    
    ```bash
    bcftools view 5NM_2Kb94aep1.vcf  -Oz -o 5NM_2Kb94aep1.vcf.gz
    ```
    ```
    bcftools view 5NM_2Kb94_k35aep1.vcf -Oz -o 5NM_2Kb94_k35aep1.vcf.gz
    ```
    ```
    bcftools index 5NM_2Kb94aep1.vcf.gz
    ```
    ``` 
    bcftools index 5NM_2Kb94_k35aep1.vcf.gz
    ```
    ```
    bcftools isec 5NM_2Kb94aep1.vcf.gz 5NM_2Kb94_k35aep1.vcf.gz -p isec_5NM_2Kb94diff_k
    ```

## Extract distance among paths

!!! terminal "code"

    ```bash

    odgi paths -i 5NM_2Kb94.gfa -d -D 'AAAA' >5NM_2Kb94.gfa_similarity
    cut -f 1,2,6 5NM_2Kb94.gfa_similarity>5NM_2Kb94.gfa_similarity_cut

    ```
!!! terminal "code"

    ```bash
  
    #Using R for distanc clustering
    module purge
    module load R/4.0.1-gimkl-2020a
    R
    ```
!!! r-project "code"
    
    ```r linenums="1"
    library(reshape)
    library(ape)
 
    # read in the data
    dat=read.csv("./5NM_2Kb94.gfa_similarity_cut",sep="\t")
    dat
    # use reshape's cast function to change to matrix
    m <- cast(dat, group.a ~ group.b)
    m
    # set the row names
    rownames(m) <- m[,1]
    rownames(m)
 
 
    # change the matrix to a distance matrix
    d <- dist(m)
    d
 
    # do hierarchical clustering
    h <- hclust(d)
 
    h
    # plot the dendrogram
    #plot(h)
 
    # use ape's as phylo function
    tree <- as.phylo(h)
    # export as newick for viewing in figtree
    write.tree(phy=tree, file = '5NM_2Kb94_distance.tree')
    ```


