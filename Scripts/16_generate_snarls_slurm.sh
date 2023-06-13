#!/bin.bash

#SBATCH --account       ga03793
#SBATCH --job-name      vgmap_generate_snarls
#SBATCH --cpus-per-task 8
#SBATCH --mem           4G
#SBATCH --time          1:00:00

module purge
module load vg/1.46.0


input_folder=/home/zyang/pg_workshop/graph_NGS/refs
output_folder=/home/zyang/pg_workshop/graph_NGS/refs

vg snarls $input_folder/4Sim_1K96_256.xg > $output_folder/4Sim_1K96_256.xg.snarls
