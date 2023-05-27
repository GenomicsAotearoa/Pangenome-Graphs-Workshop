### A script to simulate and compare variant calls generated using linear method (`bwa mem`) and graph method (`vg giraffe`)

This script peforms the below tasks; (_Nesi folder : /nesi/nobackup/ga03793/pg_workshop/vc_compare_script/_)
1. Simulate new sequence with predefined SNPs and INDELs count usnig [simuG](https://github.com/yjx1217/simuG) and create ground truth VCF file (call it as `groud_truth.vcf`)
2. Simulate reads from the new sequence with with specific coverage depth and read length using `wgsim`
3. Map the reads with the reference using `bwa mem` and generate VCF file (call it as `simulated.vcf`)
4. Compare groud_truth.vcf and simulated.vcf using `bcftools isec` and generate a report
5. Repeat steps 3 and 4 using `vg giraffe`
6. Generate a comparision stats report

The script [sim_vc_compare.sh](https://github.com/nuzla/Pangenome-Graphs-Workshop/blob/main/Scripts/sim_vc_compare.sh) accept following options. 

````
$ ./sim_vc_compare.sh -h
Program : sim_vc_compare
Version : 1.0
Usage   : sim_vc_compare.sh [options]
Options :
-r | --ref STR reference sequence file
-s | --snp INT Number of SNPs to simulate
-i | --indel INT Number of INDELs to simulate
-d | --depth INT Cover depth of the reads
-l | --length INT of a read
-o | --output STR Output folder name
-h | --help Display this help message
````

For an example if we try the below 
```
./sim_vc_compare.sh --ref GCF_000191525.1_ASM19152v1_genomic.fna -snp 5000 --indel 1000
```
Script will produce the below report for `bwa mem`.

```
+------------------------------------------------+
|  REPORT (BWA MEM)                              |
+------------------------------------------------+
|  Ground Truth SNPs                =      5,000 |
|  Ground Truth INDELs              =      1,000 |
|  Identified SNPs in Simulation    =      6,646 |
|  Identified INDELs in Simulation  =      1,666 |
|  SNPs Private to Simulation       =      1,831 |
|  INDELs Private to Simulation     =      1,153 |
|  Exact Matched SNPs               =      4,815 |
|  Exact Matched INDELs             =        513 |
|  True Positive (TP)               =      5,328 |
|  False Positive (FP)              =      2,984 |
|  True Negative (TN)               =  2,239,982 |
|  False Negative (FN)              =        672 |
+------------------------------------------------+
|  Sensitivity                      =   88.8000% |
|  Specificity                      =   99.8669% |
|  F1 Score                         =   74.4550% |
+------------------------------------------------+
```

#### Definitions
1. True Positive (TP) = SNPs+INDELs which are exactly matched in `groud_truth.vcf` and `simulated.vcf` (TP=4,815+513=5,328)
2. False Positive (FP) = SNPs+INDELs private to `simulated.vcf` and not found in `groud_truth.vcf` (FP=1,831+1,153=2,984)
3. True Negative (TN) = Length of Reference the Sequence - Ground Truth SNPs - Ground Truth INDELs - False Positive. (TN=2,248,966-5,000-1,000-2,984=2,239,982)
4. False Negative (FN) = Ground Truth SNPs+Ground Truth INDELs - True Positive. (FN=5,000+1,000-5,328=672)
5. Sensitivity and Specificity will be, 

```math
\begin{aligned}
Sensitivity  & = \frac{TP}{TP+FN} \\
              &  = \frac{5,328}{5,328+672} \\
              & = 88.8000\% \\ \\
Specificity & = \frac{TN}{TN+FP} \\
            &  = \frac{2,239,982}{2,239,982+2,984} \\
            & = 99.8669\% \\
F1\:Score & = \frac{TP}{TP+\frac{1}{2}(FP+FN)} \\
            &  = \frac{5,328}{5,328+0.5\times(2984+672)} \\ \\
            & = 74.4550\% \\
\end{aligned}
```

Script will also produce the below report for `vg giraffe` and the final comparison report.

```
+------------------------------------------------+
|  REPORT (VG GIRAFFE)                           |
+------------------------------------------------+
|  Ground Truth SNPs                =      5,000 |
|  Ground Truth INDELs              =      1,000 |
|  Identified SNPs in Simulation    =      6,832 |
|  Identified INDELs in Simulation  =      1,486 |
|  SNPs Private to Simulation       =      1,839 |
|  INDELs Private to Simulation     =        753 |
|  Exact Matched SNPs               =      4,993 |
|  Exact Matched INDELs             =        733 |
|  True Positive (TP)               =      5,726 |
|  False Positive (FP)              =      2,592 |
|  True Negative (TN)               =  2,240,374 |
|  False Negative (FN)              =        274 |
+------------------------------------------------+
|  Sensitivity                      =   95.4333% |
|  Specificity                      =   99.8844% |
|  F1 Score                         =   79.9832% |
+------------------------------------------------+
```

```
+-------------------------------------------------------------------------------------------------------------------------+
|  Method        |     TP       |     TN       |     FP       |     FN       |  Sensitivity |  Specificity |   F1 Score   |
+-------------------------------------------------------------------------------------------------------------------------+
|  bwa mem       |       5,328  |   2,239,982  |       2,984  |         672  |    88.8000%  |    99.8669%  |    74.4550%  |
+-------------------------------------------------------------------------------------------------------------------------+
|  vg giraffe    |       5,726  |   2,240,374  |       2,592  |         274  |    95.4333%  |    99.8844%  |    79.9832%  |
+-------------------------------------------------------------------------------------------------------------------------+
```

#### _Reference_
[Alignment of High-Throughput Sequencing Data using BWA](https://ucdavis-bioinformatics-training.github.io/2017-August-Variant-Analysis-Workshop/wednesday/alignment.html)
[Mapping short reads with Giraffe](https://github.com/vgteam/vg/wiki/Mapping-short-reads-with-Giraffe)


