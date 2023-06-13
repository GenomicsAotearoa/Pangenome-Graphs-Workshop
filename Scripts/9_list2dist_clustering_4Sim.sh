#!/bin.bash

#SBATCH --account       ga03793
#SBATCH --job-name      4Sim_1K96_distance_clustering
#SBATCH --cpus-per-task 8 
#SBATCH --mem           4G
#SBATCH --time          1:00:00

module purge
module load R/4.0.1-gimkl-2020a

Rscript 8_list2dist_clustering_4Sim.R
