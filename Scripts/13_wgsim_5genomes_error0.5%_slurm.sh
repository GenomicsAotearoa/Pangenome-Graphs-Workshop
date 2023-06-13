#!/bin.bash

#SBATCH --account       ga03793
#SBATCH --job-name      5NGS_simulate
#SBATCH --cpus-per-task 8
#SBATCH --mem           4G
#SBATCH --time          1:00:00

module purge
module load SAMtools/1.16.1-GCC-11.3.0
 
input_folder=/home/zyang/pg_workshop/dataset_for_pg_workshop/12_genomes_for_NGS_simulation
output_folder=/home/zyang/pg_workshop/graph_NGS/simu_NGS_data

for f in NC_017518_6k.fa ST154_6k.fa ST154Sim_6k.fa ST41Sim_6k.fa ST42Sim_6k.fa 

do

x=$(basename $f .fa)
echo ${x}

wgsim $input_folder/${x}.fa -N 1000000 -1 150 -2 150  -e 0.005 -r 0 -R 0 -X 0 $output_folder/${x}.wgsim_er0.005.R1.fq  $output_folder/${x}.wgsim_er0.005.R2.fq
gzip $output_folder/${x}.wgsim_er0.005.R1.fq
gzip $output_folder/${x}.wgsim_er0.005.R2.fq

done
