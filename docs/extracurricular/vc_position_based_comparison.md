# How to compare 2 VCF files based on the exact variant position?
When simulating the new sequences in the previous procedure [Preparing required data files with simuG tool](./preparing_data_files.md) related VCF files were also generated (e.g. `Simulation_SNP_5000.refseq2simseq.SNP.vcf`). We can compare that VCF file with linear based VCF file and the graph based files to find variant differences considering the exact bp positions. `bcftools isec` command can be used for this. 

If the ground truth simulated VCF is `Simulation_SNP_5000.refseq2simseq.SNP.vcf` (generated [here](./preparing_data_files.md)) and the linear method based vcf file (generated [here](./mapping_reads.md)) is `Simulation_SNP_5000.bwa.30x.100R.vcf`; first we need to `bgzip` the files and make indexes. 

### Ground Truth vs Linear based variant comparison 
<!--This seems redundant, the earlier exercise already generated the tables-->
!!! terminal "code"
    
    ```bash
    bgzip Simulation_SNP_5000.refseq2simseq.SNP.vcf
    bgzip Simulation_SNP_5000.bwa.30x.100R.vcf
    bcftools index Simulation_SNP_5000.bwa.30x.100R.vcf.gz
    bcftools index Simulation_SNP_5000.refseq2simseq.SNP.vcf.gz
    ```

Then use the below commands to compare. Here the option `-c none` will do the exact matching `-p` for specifying the output folder. 

!!! terminal "code"

    ```bash
    bcftools isec -c none -p output_SNP_5000_bwa.30x.100R Simulation_SNP_5000.refseq2simseq.SNP.vcf.gz Simulation_SNP_5000.bwa.30x.100R.vcf.gz
    
    ls -1sh output_SNP_5000_bwa.30x.100R
    ```

!!! success "Output"

    ```
    total 2.3M
    256K 0000.vcf
    512K 0001.vcf
    768K 0002.vcf
    768K 0003.vcf
     512 README.txt
    ```

The `README.txt` file has the details about what each vcf file has. 

!!! terminal "code"

    ```bash
    cat output_SNP_5000_bwa.30x.100R/README.txt
    ```

!!! success "Output"

    ```
    This file was produced by vcfisec.
    The command line was:   bcftools isec  -c none -p output_SNP_5000_bwa.30x.100R Simulation_SNP_5000.refseq2simseq.SNP.vcf.gz Simulation_SNP_5000.bwa.30x.100R.vcf.gz
    
    Using the following file names:
    output_SNP_5000_bwa.30x.100R/0000.vcf   for records private to  Simulation_SNP_5000.refseq2simseq.SNP.vcf.gz
    output_SNP_5000_bwa.30x.100R/0001.vcf   for records private to  Simulation_SNP_5000.bwa.30x.100R.vcf.gz
    output_SNP_5000_bwa.30x.100R/0002.vcf   for records from Simulation_SNP_5000.refseq2simseq.SNP.vcf.gz shared by both    Simulation_SNP_5000.refseq2simseq.SNP.vcf.gz Simulation_SNP_5000.bwa.30x.100R.vcf.gz
    output_SNP_5000_bwa.30x.100R/0003.vcf   for records from Simulation_SNP_5000.bwa.30x.100R.vcf.gz shared by both Simulation_SNP_5000.refseq2simseq.SNP.vcf.gz Simulation_SNP_5000.bwa.30x.100R.vcf.gz
    ```
    
As per this explanation, `0002.vcf` file should have the True Positive (TP) stats. 

!!! terminal "code"

    ```bash
    bcftools stats output_SNP_5000_bwa.30x.100R/0002.vcf | head -40
    ```

    ??? success "Output"
        ```bash
        # This file was produced by bcftools stats (1.9+htslib-1.9) and can be plotted using plot-vcfstats.
        # The command line was: bcftools stats  output_SNP_5000_bwa.30x.100R/0002.vcf
        #
        # Definition of sets:
        # ID    [2]id   [3]tab-separated file names
        ID      0       output_SNP_5000_bwa.30x.100R/0002.vcf
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
        SN      0       number of samples:      0
        SN      0       number of records:      4667
        SN      0       number of no-ALTs:      0
        SN      0       number of SNPs: 4667
        SN      0       number of MNPs: 0
        SN      0       number of indels:       0
        SN      0       number of others:       0
        SN      0       number of multiallelic sites:   0
        SN      0       number of multiallelic SNP sites:       0
        # TSTV, transitions/transversions:
        # TSTV  [2]id   [3]ts   [4]tv   [5]ts/tv        [6]ts (1st ALT) [7]tv (1st ALT) [8]ts/tv (1st ALT)
        TSTV    0       1550    3117    0.50    1550    3117    0.50
        # SiS, Singleton stats:
        # SiS   [2]id   [3]allele count [4]number of SNPs       [5]number of transitions        [6]number of transversions      [7]number of indels     [8]repeat-consistent    [9]repeat-inconsistent        [10]not applicable
        SiS     0       1       4667    1550    3117    0       0       0       0
        # AF, Stats by non-reference allele frequency:
        # AF    [2]id   [3]allele frequency     [4]number of SNPs       [5]number of transitions        [6]number of transversions      [7]number of indels     [8]repeat-consistent [9]repeat-inconsistent   [10]not applicable
        AF      0       0.000000        4667    1550    3117    0       0       0       0
        ```
    
When simulating multiple types of variants using `simuG` tool it creates multiple VCF files as well. We need to merge them and create one VCF file before comparing. For an example for the simulated sample NC_017518.1_SNP_4000_INDEL_4000, 

!!! terminal "code"

    ```bash
    ls -1hs Simulation_SNP_4000_INDEL_4000.*
    ```

!!! success "Output"

    ```
    768K Simulation_SNP_4000_INDEL_4000.refseq2simseq.INDEL.vcf
    768K Simulation_SNP_4000_INDEL_4000.refseq2simseq.map.txt
    768K Simulation_SNP_4000_INDEL_4000.refseq2simseq.SNP.vcf
    2.3M Simulation_SNP_4000_INDEL_4000.simseq.genome.fa
    ```

<!--This seems redundant, the earlier exercise already generated the tables-->
zip the files and making indexes. 

!!! terminal "code"
    ```bash
    bgzip Simulation_SNP_4000_INDEL_4000.refseq2simseq.INDEL.vcf
    bgzip Simulation_SNP_4000_INDEL_4000.refseq2simseq.SNP.vcf
    bcftools index Simulation_SNP_4000_INDEL_4000.refseq2simseq.SNP.vcf.gz
    bcftools index Simulation_SNP_4000_INDEL_4000.refseq2simseq.INDEL.vcf.gz
    ```

Merge 2 vcf files using `bcftools`

!!! terminal "code"
    ```bash
    bcftools merge Simulation_SNP_4000_INDEL_4000.refseq2simseq.SNP.vcf.gz Simulation_SNP_4000_INDEL_4000.refseq2simseq.INDEL.vcf.gz -O z -o Simulation_SNP_4000_INDEL_4000.vcf.gz
    ```

Stats of the merged vcf file should show 4000 SNPs and 4000 INDELs.

!!! terminal "code"
    ```bash
    bcftools stats Simulation_SNP_4000_INDEL_4000.vcf.gz | less
    ```
    ??? success "Output"
        ```
        # This file was produced by bcftools stats (1.9+htslib-1.9) and can be plotted using plot-vcfstats.
        # The command line was: bcftools stats  Simulation_SNP_4000_INDEL_4000.vcf.gz
        #
        # Definition of sets:
        # ID    [2]id   [3]tab-separated file names
        ID      0       Simulation_SNP_4000_INDEL_4000.vcf.gz
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
        SN      0       number of samples:      0
        SN      0       number of records:      8000
        SN      0       number of no-ALTs:      0
        SN      0       number of SNPs: 4000
        SN      0       number of MNPs: 0
        SN      0       number of indels:       4000
        SN      0       number of others:       0
        SN      0       number of multiallelic sites:   0
        SN      0       number of multiallelic SNP sites:       0
        # TSTV, transitions/transversions:
        # TSTV  [2]id   [3]ts   [4]tv   [5]ts/tv        [6]ts (1st ALT) [7]tv (1st ALT) [8]ts/tv (1st ALT)
        TSTV    0       1372    2628    0.52    1372    2628    0.52
        # SiS, Singleton stats:
        # SiS   [2]id   [3]allele count [4]number of SNPs       [5]number of transitions        [6]number of transversions      [7]number of indels     [8]repeat-consistent    [9]repeat-inconsistent  [10]not applicable
        SiS     0       1       4000    1372    2628    4000    0       0       4000
        ```

Now we can apply the same procedure for comparison using `bcftools isec`. <!-- This suggests some iteration, manual or looped, but I'm not sure where the iteration is implemented here? -->

<!-- 
What are learners looking for here? What is the downstream context here?
-->

<!---

### Ground Truth vs Graph based variant comparision 
As explain in the main page we have applied some filtering to exatract the PGGB VCF and then compare to the Ground Truth VCF using bcftool isec command.

```bash
#for 95
bcftools view -s NC_017518.1_SNP_5000 --min-ac=1 pggb_1k95_NC_017518.1.vcf  >  Simulation_SNP_5000.graph95.vcf 
bgzip Simulation_SNP_5000.graph95.vcf
bcftools index Simulation_SNP_5000.graph95.vcf.gz 
bcftools isec -c none -p output_SNP_5000_graph95 Simulation_SNP_5000.refseq2simseq.SNP.vcf.gz Simulation_SNP_5000.graph95.vcf.gz

$ cd output_SNP_5000_graph95

$ ls -1hs 
total 5.1M
256K 0000.vcf
3.3M 0001.vcf
768K 0002.vcf
768K 0003.vcf
 512 README.txt
 
$ cat README.txt 
This file was produced by vcfisec.
The command line was:   bcftools isec  -c none -p output_SNP_5000_graph95 Simulation_SNP_5000.refseq2simseq.SNP.vcf.gz Simulation_SNP_5000.graph95.vcf.gz

Using the following file names:
output_SNP_5000_graph95/0000.vcf        for records private to  Simulation_SNP_5000.refseq2simseq.SNP.vcf.gz
output_SNP_5000_graph95/0001.vcf        for records private to  Simulation_SNP_5000.graph95.vcf.gz
output_SNP_5000_graph95/0002.vcf        for records from Simulation_SNP_5000.refseq2simseq.SNP.vcf.gz shared by both    Simulation_SNP_5000.refseq2simseq.SNP.vcf.gz Simulation_SNP_5000.graph95.vcf.gz
output_SNP_5000_graph95/0003.vcf        for records from Simulation_SNP_5000.graph95.vcf.gz shared by both      Simulation_SNP_5000.refseq2simseq.SNP.vcf.gz Simulation_SNP_5000.graph95.vcf.gz

$ bcftools stats 0002.vcf 
# This file was produced by bcftools stats (1.9+htslib-1.9) and can be plotted using plot-vcfstats.
# The command line was: bcftools stats  0002.vcf
#
# Definition of sets:
# ID    [2]id   [3]tab-separated file names
ID      0       0002.vcf
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
SN      0       number of samples:      0
SN      0       number of records:      4838
SN      0       number of no-ALTs:      0
SN      0       number of SNPs: 4838
SN      0       number of MNPs: 0
SN      0       number of indels:       0
SN      0       number of others:       0
SN      0       number of multiallelic sites:   0
SN      0       number of multiallelic SNP sites:       0
# TSTV, transitions/transversions:
# TSTV  [2]id   [3]ts   [4]tv   [5]ts/tv        [6]ts (1st ALT) [7]tv (1st ALT) [8]ts/tv (1st ALT)
TSTV    0       1601    3237    0.49    1601    3237    0.49
# SiS, Singleton stats:
# SiS   [2]id   [3]allele count [4]number of SNPs       [5]number of transitions        [6]number of transversions      [7]number of indels     [8]repeat-consistent    [9]repeat-inconsistent [10]not applicable
SiS     0       1       4838    1601    3237    0       0       0       0
# AF, Stats by non-reference allele frequency:
# AF    [2]id   [3]allele frequency     [4]number of SNPs       [5]number of transitions        [6]number of transversions      [7]number of indels     [8]repeat-consistent    [9]repeat-inconsistent [10]not applicable
AF      0       0.000000        4838    1601    3237    0       0       0       0
# QUAL, Stats by quality:
# QUAL  [2]id   [3]Quality      [4]number of SNPs       [5]number of transitions (1st ALT)      [6]number of transversions (1st ALT)    [7]number of indels
QUAL    0       998     4838    1601    3237    0
# IDD, InDel distribution:
# IDD   [2]id   [3]length (deletions negative)  [4]count
# ST, Substitution types:
# ST    [2]id   [3]type [4]count
ST      0       A>C     391
ST      0       A>G     406
ST      0       A>T     420
ST      0       C>A     410
ST      0       C>G     429
ST      0       C>T     396
ST      0       G>A     380
ST      0       G>C     409
ST      0       G>T     391
ST      0       T>A     392
ST      0       T>C     419
ST      0       T>G     395
# DP, Depth distribution
# DP    [2]id   [3]bin  [4]number of genotypes  [5]fraction of genotypes (%)    [6]number of sites      [7]fraction of sites (%)
```
Using the same procedure to follow PGGB 98 files.
```
#for 98
bcftools view -s NC_017518.1_SNP_5000 --min-ac=1 pggb_1k98_NC_017518.1.vcf  >  Simulation_SNP_5000.graph98.vcf
bgzip Simulation_SNP_5000.graph98.vcf
bcftools index Simulation_SNP_5000.graph98.vcf.gz 
bcftools isec -c none -p output_SNP_5000_graph98 Simulation_SNP_5000.refseq2simseq.SNP.vcf.gz Simulation_SNP_5000.graph98.vcf.gz

$ cd output_SNP_5000_graph98
$ ls -1hs
total 5.1M
256K 0000.vcf
3.3M 0001.vcf
768K 0002.vcf
768K 0003.vcf
 512 README.txt

$ cat README.txt 
This file was produced by vcfisec.
The command line was:   bcftools isec  -c none -p output_SNP_5000_graph98 Simulation_SNP_5000.refseq2simseq.SNP.vcf.gz Simulation_SNP_5000.graph98.vcf.gz

Using the following file names:
output_SNP_5000_graph98/0000.vcf        for records private to  Simulation_SNP_5000.refseq2simseq.SNP.vcf.gz
output_SNP_5000_graph98/0001.vcf        for records private to  Simulation_SNP_5000.graph98.vcf.gz
output_SNP_5000_graph98/0002.vcf        for records from Simulation_SNP_5000.refseq2simseq.SNP.vcf.gz shared by both    Simulation_SNP_5000.refseq2simseq.SNP.vcf.gz Simulation_SNP_5000.graph98.vcf.gz
output_SNP_5000_graph98/0003.vcf        for records from Simulation_SNP_5000.graph98.vcf.gz shared by both      Simulation_SNP_5000.refseq2simseq.SNP.vcf.gz Simulation_SNP_5000.graph98.vcf.gz

$ bcftools stats 0002.vcf 
# This file was produced by bcftools stats (1.9+htslib-1.9) and can be plotted using plot-vcfstats.
# The command line was: bcftools stats  0002.vcf
#
# Definition of sets:
# ID    [2]id   [3]tab-separated file names
ID      0       0002.vcf
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
SN      0       number of samples:      0
SN      0       number of records:      4843
SN      0       number of no-ALTs:      0
SN      0       number of SNPs: 4843
SN      0       number of MNPs: 0
SN      0       number of indels:       0
SN      0       number of others:       0
SN      0       number of multiallelic sites:   0
SN      0       number of multiallelic SNP sites:       0
# TSTV, transitions/transversions:
# TSTV  [2]id   [3]ts   [4]tv   [5]ts/tv        [6]ts (1st ALT) [7]tv (1st ALT) [8]ts/tv (1st ALT)
TSTV    0       1602    3241    0.49    1602    3241    0.49
# SiS, Singleton stats:
# SiS   [2]id   [3]allele count [4]number of SNPs       [5]number of transitions        [6]number of transversions      [7]number of indels     [8]repeat-consistent    [9]repeat-inconsistent  [10]not applicable
SiS     0       1       4843    1602    3241    0       0       0       0
# AF, Stats by non-reference allele frequency:
# AF    [2]id   [3]allele frequency     [4]number of SNPs       [5]number of transitions        [6]number of transversions      [7]number of indels     [8]repeat-consistent    [9]repeat-inconsistent  [10]not applicable
AF      0       0.000000        4843    1602    3241    0       0       0       0
# QUAL, Stats by quality:
# QUAL  [2]id   [3]Quality      [4]number of SNPs       [5]number of transitions (1st ALT)      [6]number of transversions (1st ALT)    [7]number of indels
QUAL    0       998     4843    1602    3241    0
# IDD, InDel distribution:
# IDD   [2]id   [3]length (deletions negative)  [4]count
# ST, Substitution types:
# ST    [2]id   [3]type [4]count
ST      0       A>C     392
ST      0       A>G     406
ST      0       A>T     419
ST      0       C>A     412
ST      0       C>G     429
ST      0       C>T     398
ST      0       G>A     381
ST      0       G>C     411
ST      0       G>T     391
ST      0       T>A     392
ST      0       T>C     417
ST      0       T>G     395
# DP, Depth distribution
# DP    [2]id   [3]bin  [4]number of genotypes  [5]fraction of genotypes (%)    [6]number of sites      [7]fraction of sites (%)
```
--->
