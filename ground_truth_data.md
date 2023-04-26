# Script for finding Ground Truth SNP count

```bash
#Load required modules with specific versions
module purge
module load BCFtools/1.9-GCC-7.4.0
module load SAMtools/1.9-GCC-7.4.0
module load BWA/0.7.17-GCC-9.2.0

mkdir SNP

#Splitting the 4Sim.fa dataset into four files
samtools faidx 4Sim.fa NC_neisseria > ./SNP/NC_neisseria.fa
samtools faidx 4Sim.fa Sim1_3k > ./SNP/Sim1_3k.fa
samtools faidx 4Sim.fa Sim2_4k > ./SNP/Sim2_4k.fa
samtools faidx 4Sim.fa Sim3_5k > ./SNP/Sim3_5k.fa
cd SNP/

#indexing the reference sample
bwa index NC_neisseria.fa 

#creating VCF files for each sample file considering NC_neisseria.fa as reference
bwa mem -R "@RG\tID:Sim1_3k\tSM:Sim1_3k\tLB:L1" NC_neisseria.fa Sim1_3k.fa > Sim1_3k.sam
samtools view -bS Sim1_3k.sam | samtools sort - > Sim1_3k.bam
bcftools mpileup -Ou -f NC_neisseria.fa Sim1_3k.bam | bcftools call -vmO z -o Sim1_3k.vcf.gz
bcftools index Sim1_3k.vcf.gz 

bwa mem -R "@RG\tID:Sim2_4k\tSM:Sim2_4k\tLB:L1" NC_neisseria.fa Sim2_4k.fa > Sim2_4k.sam
samtools view -bS Sim2_4k.sam | samtools sort - > Sim2_4k.bam
bcftools mpileup -Ou -f NC_neisseria.fa Sim2_4k.bam | bcftools call -vmO z -o Sim2_4k.vcf.gz
bcftools index Sim2_4k.vcf.gz 

bwa mem -R "@RG\tID:Sim3_5k\tSM:Sim3_5k\tLB:L1" NC_neisseria.fa Sim3_5k.fa > Sim3_5k.sam
samtools view -bS Sim3_5k.sam | samtools sort - > Sim3_5k.bam
bcftools mpileup -Ou -f NC_neisseria.fa Sim3_5k.bam | bcftools call -vmO z -o Sim3_5k.vcf.gz
bcftools index Sim3_5k.vcf.gz 

#merging the sample three VCF outputs one VCF file
bcftools merge -O v -o 4Sim.vcf Sim1_3k.vcf.gz Sim2_4k.vcf.gz Sim3_5k.vcf.gz

#Find the variants stats
bcftools stats 4Sim.vcf > 4Sim.vcf.stats

```
