#!/usr/bin/bash

#SBATCH --account       ga03793
#SBATCH --job-name      4Sim_1K96_K79
#SBATCH --cpus-per-task 8 
#SBATCH --mem           4G
#SBATCH --time          1:00:00

module purge
module load Singularity

#export container to a variable for convenience
WD=/home/zyang/pg_workshop #Working Directory
container=/nesi/project/ga03793/software/pggb/pggb_0.5.3.simg
data=/home/zyang/pg_workshop/4Sim.fa
output=/home/zyang/pg_workshop


#Bind filesystem to container image 
export SINGULARITY_BIND="${WD}, /nesi/project/ga03793/"

singularity exec ${container} pggb -i $data -s 1000 -p 96 -n 4 -K 79 -t 24 -S -m -o $output/4Sim_1K96_K79 -V 'NC_017518:#' 
