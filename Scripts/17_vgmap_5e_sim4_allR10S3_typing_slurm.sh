#!/bin.bash

#SBATCH --account       ga03793
#SBATCH --job-name      5e_vgmap_genotying 
#SBATCH --cpus-per-task 8
#SBATCH --mem           4G
#SBATCH --time          24:00:00

module purge
module load vg/1.46.0

data_gam=/home/zyang/pg_workshop/graph_NGS/graph_based_mapping/*.wgsim_er0.005.vgmap_4Sim.gam
input=/home/zyang/pg_workshop/graph_NGS/graph_based_mapping
output=/home/zyang/pg_workshop/graph_NGS/vgmap_12e_sim4_allR10S3_typing
graph_xg=/home/zyang/pg_workshop/graph_NGS/refs/4Sim_1K96_256.xg
snarls_file=/home/zyang/pg_workshop/graph_NGS/refs/4Sim_1K96_256.xg.snarls



#compute snarls
#vg snarls $graph_xg >$snarls_file

for f in $data_gam
do 

x=$(basename $f .wgsim_er0.005.vgmap_4Sim.gam)
echo ${x}


#Calculate the surpport reads ingoring mapping and base quality <5
#vg pack -t 48 -x $graph_xg -g $input/${x}.wgsim_er0.005.vgmap_4Sim.gam -Q 5 -o $output/${x}vgmap_Sim4_256_aln.pack

#Calculate the surpport reads
vg pack -t 12 -x $graph_xg -g $input/${x}.wgsim_er0.005.vgmap_4Sim.gam -o $output/${x}vgmap_sim4_256_aln.pack

#call variant using the same coordinates and including reference calls (for following compare)
vg call -t 12 -m 3,10 $graph_xg -k $output/${x}vgmap_sim4_256_aln.pack -r $snarls_file -a  >$output/${x}vgmap_sim4_256_aln.pack_allR10S3.vcf 

done
