#!/bin/bash

#SBATCH --account       ga03793
#SBATCH --job-name      4Sim_1K96_K79
#SBATCH --cpus-per-task 8 
#SBATCH --mem           4G
#SBATCH --time          1:00:00

module purge



WD=/home/zyang/pg_workshop #Working Directory

data=/home/zyang/pg_workshop/4Sim.fa
output=/home/zyang/pg_workshop


 


singularity exec ${container} pggb -i $data -s 1000 -p 96 -n 4 -K 79 -t 24 -S -m -o $output/4Sim_1K96_K79 -V 'NC_017518:#' 
