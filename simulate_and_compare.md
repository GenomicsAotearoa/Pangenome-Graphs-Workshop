### A script to simulate and compare variant calls generated using linear method (`bwa mem`)

This script peforms the below tasks; (_Nesi folder : /nesi/nobackup/ga03793/pg_workshop/vc_compare_script/_)
1. Simulate new sequence with predefined SNPs and INDELs count usnig [simuG](https://github.com/yjx1217/simuG) and create ground truth VCF file (call it as `groud_truth.vcf`)
2. Simulate reads from the new sequence with with specific coverage depth and read length using `wgsim`
3. Map the reads with the reference using `bwa mem` and generate VCF file (call it as `simulated.vcf`)
4. Compare groud_truth.vcf and simulated.vcf using `bcftools isec` and generate a report

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
Script will produce the below final report.

```
--------------------------------------------------
|  REPORT                                        |
--------------------------------------------------
|  Ground Truth SNPs                =      5,000 |
|  Ground Truth INDELs              =      1,000 |
|  Identified SNPs in Simulation    =      6,640 |
|  Identified INDELs in Simulation  =      1,690 |
|  SNPs Private to Simulation       =      1,821 |
|  INDELs Private to Simulation     =      1,176 |
|  Exact Matched SNPs               =      4,819 |
|  Exact Matched INDELs             =        514 |
|  True Positive (TP)               =      5,333 |
|  False Positive (FP)              =      2,997 |
|  True Negative (TN)               =  2,239,969 |
|  False Negative (FN)              =        667 |
|------------------------------------------------
|  Sensitivity                      =   88.8833% |
|  Specificity                      =   99.8663% |
--------------------------------------------------
```

#### Definitions
1. True Positive (TP) = SNPs+INDELs which are exactly matched in `groud_truth.vcf` and `simulated.vcf` (TP=4,819+514=5,333)
2. False Positive (FP) = SNPs+INDELs private to `simulated.vcf` and not found in `groud_truth.vcf` (FP=1,821+1,176=2,997)
3. True Negative (TN) = Length of Reference the Sequence - Ground Truth SNPs - Ground Truth INDELs - False Positive. (TN=2,248,966-5,000-1,000-2,997=2,239,969)
4. False Negative (FN) = Ground Truth SNPs+Ground Truth INDELs - True Positive. (FN=5,000+1,000-5,333=667)
5. Sensitivity and Specificity will be, 

```math
\begin{aligned}
Sensitivity  & = \frac{TP}{TP+FN} \\
              &  = \frac{5,333}{5,333+667} \\
              & = 88.8833\% \\ \\
Specificity & = \frac{TN}{TN+FP} \\
            &  = \frac{2,239,969}{2,239,969+2,997} \\
            & = 99.8663\% \\
\end{aligned}
```

