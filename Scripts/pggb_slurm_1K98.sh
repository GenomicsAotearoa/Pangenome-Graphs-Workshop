#!/bin/bash

#SBATCH --account       ga03793
#SBATCH --job-name      NC_017518.1_1K98
#SBATCH --cpus-per-task 8 
#SBATCH --mem           4G
#SBATCH --time          1:00:00

module purge



WD=/nesi/nobackup/nesi02659/pg_workshop #Working Directory

data=${WD}/ASM19152v1_pgsim.fa

 


singularity exec ${container} pggb -i $data -s 1000 -p 98 -n 6 -k 79 -t 24 -S -m -o output_1K98 -V 'NC_017518.1:#' 
