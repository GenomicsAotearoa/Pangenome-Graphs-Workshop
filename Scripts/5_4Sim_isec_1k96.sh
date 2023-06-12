#!/usr/bin/bash

#SBATCH --account       ga03793
#SBATCH --job-name      4Sim_1k96_isec
#SBATCH --cpus-per-task 8
#SBATCH --mem           4G
#SBATCH --time          1:00:00

module purge
module load BCFtools/1.16-GCC-11.3.0

input_folder=/home/zyang/pg_workshop/vg_deconstruct
output_folder=/home/zyang/pg_workshop/vg_deconstruct


bcftools view $input_folder/4Sim_1K96aep1.vcf  -Oz -o $output_folder/4Sim_1K96aep1.vcf.gz
bcftools view $input_folder/4Sim_1K96_K79aep1.vcf -Oz -o $output_folder/4Sim_1K96_K79aep1.vcf.gz

bcftools index $output_folder/4Sim_1K96aep1.vcf.gz
bcftools index $output_folder/4Sim_1K96_K79aep1.vcf.gz

bcftools isec $output_folder/4Sim_1K96aep1.vcf.gz $output_folder/4Sim_1K96_K79aep1.vcf.gz -p $output_folder/isec_4Sim_1K96

