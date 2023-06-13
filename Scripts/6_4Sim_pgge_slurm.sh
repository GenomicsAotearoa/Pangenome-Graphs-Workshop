#!/bin.bash

#SBATCH --account       ga03793
#SBATCH --job-name      4Sim_pgge
#SBATCH --cpus-per-task 8
#SBATCH --mem           4G
#SBATCH --time          1:00:00

module purge
module load Singularity

#export container to a variable for convenience

container=/nesi/project/nesi02659/software/pgge/pgge_032023.simg


WD=/home/zyang/pg_workshop/4Sim_pgge #Working Directory
inputGFA=/home/zyang/pg_workshop/4Sim_pgge/*.gfa
input_folder=/home/zyang/pg_workshop/4Sim_pgge
inputfa=/home/zyang/pg_workshop/4Sim.fa
output=/home/zyang/pg_workshop/4Sim_pgge
beehave=/home/zyang/pg_workshop/beehave.R


#Bind filesystem to container image
export SINGULARITY_BIND="${WD}, /nesi/project/nesi02659/"


for x in $inputGFA

do
singularity exec ${container} pgge -g $x -f $inputfa -o $output -r $beehave -b $output/pgge_4Sim_peanut_bed -l 100000 -s 5000 -t 16

done
