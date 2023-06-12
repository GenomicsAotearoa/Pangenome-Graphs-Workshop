#!/usr/bin/bash

#SBATCH --account       ga03793
#SBATCH --job-name      build_index_for_4SimGraph
#SBATCH --cpus-per-task 8
#SBATCH --mem           4G
#SBATCH --time          1:00:00

module purge
module load vg/1.46.0

cd /home/zyang/pg_workshop/graph_NGS/refs 
data=/home/zyang/pg_workshop/graph_NGS/refs/*.gfa  
tem_dir=/home/zyang/pg_workshop/graph_NGS/refs/temp_dir


for f in $data 
do 

x=$(basename $f .gfa)


#convert the graph into 256 bp chunks, saving as vg format
vg mod -X 256 ${x}.gfa >${x}_256.vg 

#build index of xg and gcsa index
vg index -b $tem_dir -t 48 -x ${x}_256.xg -g ${x}_256.gcsa -k 16 ${x}_256.vg 


#small graph is ok without prunning, complex graph will need to prune first before generating index  

### pruning: use -M if pruning fails
#vg prune -u -m node-mapping.tmp -t 48 -k 24 ${x}_256.vg > ${x}_256_chopped.vg

#vg index ${x}_256_chopped.vg -x ${x}_256_chopped.xg
### gcsa index
#vg index -b $tem_dir -t 48  -g ${x}_256_chopped.gcsa  ${x}_256_chopped.vg
done  
