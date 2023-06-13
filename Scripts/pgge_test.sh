#!/bin/bash
module load Singularity
#export container to a variable for convenience
WD=/nesi/nobackup/nesi02659/pg_workshop #Working Directory
container=/nesi/project/nesi02659/software/pgge/pgge_032023.simg
data=${WD}/4Sim.fa

#Bind filesystem to container image 
export SINGULARITY_BIND="${WD}, /nesi/project/nesi02659/"

singularity exec ${container} pgge -g ${WD}/output/*.gfa -f $data -o pgge -r ${WD}/beehave.R -b pgge/pgge_4Sim_peanut_bed -l 100000 -s 5000 -t 16
