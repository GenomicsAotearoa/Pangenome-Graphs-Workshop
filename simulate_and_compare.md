### A script for simulate and compare variant calls using linear method

This script peform the below tasks
1. Simulate new sequence with predefined SNPs and INDELs count usnig [simuG](https://github.com/yjx1217/simuG) and create ground truth VCF file (call it as groud_truth.vcf)
2. Simulate reads from the new sequence with with specific coverage depth and read length using `wgsim`
3. Map the reads with the reference using `bwa mem` and generate VCF file (call it as simulated.vcf)
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
|  Ground Truth SNPs               =      5,000  |
|  Ground Truth INDELs             =      1,000  |
|  Simulated SNPs                  =      6,717  |
|  Simulated INDELs                =      1,597  |
|  Uniq SNPs in Simulation         =      1,871  |
|  Uniq INDELs in Simulation       =      1,098  |
|  Common SNPs in Both             =      4,846  |
|  Common INDELs in Both           =        499  |
|  True Positive (TP)              =      5,345  |
|  False Positive (FP)             =      2,969  |
|  True Negative (TN)              =  2,239,997  |
|  False Negative (FN)             =        655  |
|------------------------------------------------
|  Sensitivity                     =    89.0833% |
|  Specificity                     =    99.8676% |
--------------------------------------------------
```
