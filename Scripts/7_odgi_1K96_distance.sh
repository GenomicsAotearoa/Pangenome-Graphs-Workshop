#!/usr/bin/bash

#SBATCH --account       ga03793
#SBATCH --job-name      4Sim_1K96_distance
#SBATCH --cpus-per-task 8 
#SBATCH --mem           4G
#SBATCH --time          1:00:00

module purge
module load Singularity

#export container to a variable for convenience
container=/nesi/project/nesi02659/software/odgi/odgi_0.8.2.simg
data=/home/zyang/pg_workshop/odgi_distance/4Sim_1K96.gfa
output=/home/zyang/pg_workshop/odgi_distance



singularity exec ${container} odgi paths -i $data -d -D 'AAAA' >$output/4Sim_1K96.gfa_distance
cut -f 1,2,6 $output/4Sim_1K96.gfa_distance >$output/4Sim_1K96.gfa_distance_cut 
