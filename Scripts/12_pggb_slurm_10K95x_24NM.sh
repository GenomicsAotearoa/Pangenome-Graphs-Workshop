#!/bin.bash

#SBATCH --account       ga03793
#SBATCH --job-name      24NM_10k95x
#SBATCH --cpus-per-task 8 
#SBATCH --mem           4G
#SBATCH --time          3:00:00

module purge
module load Singularity

#export container to a variable for convenience
WD=/home/zyang/pg_workshop #Working Directory
container=/nesi/project/nesi02659/software/pggb/pggb_0.5.3.simg
data=/home/zyang/pg_workshop/24NM.fa
output=/home/zyang/pg_workshop


#Bind filesystem to container image 
export SINGULARITY_BIND="${WD}, /nesi/project/nesi02659/"

singularity exec ${container} pggb -i $data -s 10000 -p 95 -n 24 -x auto -t 24 -S -m -o $output/24NM_10K95x -V 'NC_017518:#' 