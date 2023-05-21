# Script for finding Linear Reference Based Variant Calls

```bash
#Load required modules with specific versions
module purge
module load BCFtools/1.9-GCC-7.4.0
module load SAMtools/1.9-GCC-7.4.0
module load BWA/0.7.17-GCC-9.2.0

mkdir VC_compare
cd VC_compare

#Copy the below files into the folder
$ ls -1trhs
total 14M
2.3M GCF_000191525.1_ASM19152v1_genomic.fna
2.3M Simulation_INDEL_5000.simseq.genome.fa
2.5M Simulation_SNP_4000_INDEL_4000_CNV_4.simseq.genome.fa
2.3M Simulation_SNP_4000_INDEL_4000_INV_4.simseq.genome.fa
2.3M Simulation_SNP_4000_INDEL_4000.simseq.genome.fa
2.3M Simulation_SNP_5000.simseq.genome.fa

#indexing the reference sample
bwa index GCF_000191525.1_ASM19152v1_genomic.fna 

#creating VCF files for each sample file considering GCF_000191525.1_ASM19152v1_genomic.fna as reference
#export OMP_NUM_THREADS=1
bwa mem -R "@RG\tID:Simulation_SNP_5000\tSM:Simulation_SNP_5000\tLB:L1" GCF_000191525.1_ASM19152v1_genomic.fna Simulation_SNP_5000.simseq.genome.fa > Simulation_SNP_5000.sam
samtools view -bS Simulation_SNP_5000.sam | samtools sort - > Simulation_SNP_5000.bam
bcftools mpileup -Ou -f GCF_000191525.1_ASM19152v1_genomic.fna Simulation_SNP_5000.bam | bcftools call -vmO z -o Simulation_SNP_5000.vcf.gz
bcftools index Simulation_SNP_5000.vcf.gz 

bwa mem -R "@RG\tID:Simulation_INDEL_5000\tSM:Simulation_INDEL_5000\tLB:L1" GCF_000191525.1_ASM19152v1_genomic.fna Simulation_INDEL_5000.simseq.genome.fa > Simulation_INDEL_5000.sam
samtools view -bS Simulation_INDEL_5000.sam | samtools sort - > Simulation_INDEL_5000.bam
bcftools mpileup -Ou -f GCF_000191525.1_ASM19152v1_genomic.fna Simulation_INDEL_5000.bam | bcftools call -vmO z -o Simulation_INDEL_5000.vcf.gz
bcftools index Simulation_INDEL_5000.vcf.gz 

bwa mem -R "@RG\tID:Simulation_SNP_4000_INDEL_4000\tSM:Simulation_SNP_4000_INDEL_4000\tLB:L1" GCF_000191525.1_ASM19152v1_genomic.fna Simulation_SNP_4000_INDEL_4000.simseq.genome.fa > Simulation_SNP_4000_INDEL_4000.sam
samtools view -bS Simulation_SNP_4000_INDEL_4000.sam | samtools sort - > Simulation_SNP_4000_INDEL_4000.bam
bcftools mpileup -Ou -f GCF_000191525.1_ASM19152v1_genomic.fna Simulation_SNP_4000_INDEL_4000.bam | bcftools call -vmO z -o Simulation_SNP_4000_INDEL_4000.vcf.gz
bcftools index Simulation_SNP_4000_INDEL_4000.vcf.gz

bwa mem -R "@RG\tID:Simulation_SNP_4000_INDEL_4000_INV_4\tSM:Simulation_SNP_4000_INDEL_4000_INV_4\tLB:L1" GCF_000191525.1_ASM19152v1_genomic.fna Simulation_SNP_4000_INDEL_4000_INV_4.simseq.genome.fa > Simulation_SNP_4000_INDEL_4000_INV_4.sam
samtools view -bS Simulation_SNP_4000_INDEL_4000_INV_4.sam | samtools sort - > Simulation_SNP_4000_INDEL_4000_INV_4.bam
bcftools mpileup -Ou -f GCF_000191525.1_ASM19152v1_genomic.fna Simulation_SNP_4000_INDEL_4000_INV_4.bam | bcftools call -vmO z -o Simulation_SNP_4000_INDEL_4000_INV_4.vcf.gz
bcftools index Simulation_SNP_4000_INDEL_4000_INV_4.vcf.gz

bwa mem -R "@RG\tID:Simulation_SNP_4000_INDEL_4000_CNV_4\tSM:Simulation_SNP_4000_INDEL_4000_CNV_4\tLB:L1" GCF_000191525.1_ASM19152v1_genomic.fna Simulation_SNP_4000_INDEL_4000_CNV_4.simseq.genome.fa > Simulation_SNP_4000_INDEL_4000_CNV_4.sam
samtools view -bS Simulation_SNP_4000_INDEL_4000_CNV_4.sam | samtools sort - > Simulation_SNP_4000_INDEL_4000_CNV_4.bam
bcftools mpileup -Ou -f GCF_000191525.1_ASM19152v1_genomic.fna Simulation_SNP_4000_INDEL_4000_CNV_4.bam | bcftools call -vmO z -o Simulation_SNP_4000_INDEL_4000_CNV_4.vcf.gz
bcftools index Simulation_SNP_4000_INDEL_4000_CNV_4.vcf.gz

```

Now we can find the varient call stats using `bcftools stats`. 

```bash
$ bcftools stats Simulation_SNP_5000.vcf.gz | head -30
# This file was produced by bcftools stats (1.9+htslib-1.9) and can be plotted using plot-vcfstats.
# The command line was: bcftools stats  Simulation_SNP_5000.vcf.gz
#
# Definition of sets:
# ID    [2]id   [3]tab-separated file names
ID      0       Simulation_SNP_5000.vcf.gz
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
SN      0       number of records:      5000
SN      0       number of no-ALTs:      0
SN      0       number of SNPs: 5000
SN      0       number of MNPs: 0
SN      0       number of indels:       0
SN      0       number of others:       0
SN      0       number of multiallelic sites:   0
```
