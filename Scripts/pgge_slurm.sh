#!/usr/bin/bash

#SBATCH --account       ga03793
#SBATCH --job-name      4Sim_pgge
#SBATCH --cpus-per-task 8 
#SBATCH --mem           4G
#SBATCH --time          1:00:00

module purge
module load Singularity
#export container to a variable for convenience
WD=/nesi/nobackup/nesi02659/pg_workshop #Working Directory
container=/nesi/project/nesi02659/software/pgge/pgge_032023.simg
data=${WD}/4Sim.fa

#Bind filesystem to container image 
export SINGULARITY_BIND="${WD}, /nesi/project/nesi02659/"

singularity exec ${container} pgge -g ${WD}/output/*.gfa -f $data -o pgge_output -r ${WD}/beehave.R -b pgge_output/pgge_4Sim_peanut_bed -l 100000 -s 5000 -t 8 
