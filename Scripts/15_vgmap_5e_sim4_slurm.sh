#!/bin/bash

#SBATCH --account       ga03793
#SBATCH --job-name      vgmap_5e_4Sim
#SBATCH --cpus-per-task 8
#SBATCH --mem           4G
#SBATCH --time          24:00:00

module purge
module load vg/1.46.0


data=/home/zyang/pg_workshop/graph_NGS/simu_NGS_data/*R1.fq.gz
input_folder=/home/zyang/pg_workshop/graph_NGS/simu_NGS_data
output=/home/zyang/pg_workshop/graph_NGS/graph_based_mapping


index=/home/zyang/pg_workshop/graph_NGS/refs/4Sim_1K96_256.gcsa
basename=/home/zyang/pg_workshop/graph_NGS/refs/4Sim_1K96_256


for f in $data
do


x=$(basename $f R1.fq.gz)
echo ${x}

read1=${x}R1.fq.gz
read2=$(echo $read1|sed 's/R1.fq.gz/R2.fq.gz/')

echo $read2

#map paired reads using vg map
vg map -t 20  -d $basename -g $index  -f $input_folder/$read1 -f $input_folder/$read2 -N $x  > $output/${x}vgmap_4Sim.gam
#vg stats to check the mapping statistics 
vg stats -a  $output/${x}vgmap_4Sim.gam  >$output/${x}vgmap_4Sim_stats 

done


