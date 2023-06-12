#!/usr//bin/bash

#SBATCH --account       ga03793
#SBATCH --job-name      4Sim_vg_deconstruct
#SBATCH --cpus-per-task 8
#SBATCH --mem           4G
#SBATCH --time          1:00:00

module purge
module load vg/1.46.0
module load BCFtools/1.16-GCC-11.3.0

export container to a variable for convenience
inputGFA=/home/zyang/pg_workshop/vg_deconstruct/*.gfa
input_folder=/home/zyang/pg_workshop/vg_deconstruct
output=/home/zyang/pg_workshop/vg_deconstruct


for f in $inputGFA

do

x=$(basename $f .gfa)
echo ${x}

vg deconstruct -p NC_017518  -a -e $input_folder/${x}.gfa > $output/${x}aep1.vcf
bcftools stats $output/${x}aep1.vcf >$output/${x}aep1.vcf_stats 

done


