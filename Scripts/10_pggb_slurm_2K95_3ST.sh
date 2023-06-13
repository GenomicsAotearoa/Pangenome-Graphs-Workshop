#!/bin/bash

#SBATCH --account       ga03793
#SBATCH --job-name      3ST_2K95
#SBATCH --cpus-per-task 8 
#SBATCH --mem           4G
#SBATCH --time          1:00:00

module purge
module load Singularity

#export container to a variable for convenience
container=/nesi/project/nesi02659/software/pggb/pggb_0.5.3.simg
data=/home/zyang/pg_workshop/3ST.fa
output=/home/zyang/pg_workshop


singularity exec ${container} pggb -i $data -s 2000 -p 95 -n 3 -t 24 -S -m -o $output/3ST_2K95 -V 'NC_017518:#' 
