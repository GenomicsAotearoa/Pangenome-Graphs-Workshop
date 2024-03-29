#!/bin/bash

#SBATCH --account       ga03793
#SBATCH --job-name      4Sim_pgge
#SBATCH --cpus-per-task 8
#SBATCH --mem           4G
#SBATCH --time          1:00:00

module purge







WD=/home/zyang/pg_workshop/4Sim_pgge #Working Directory
inputGFA=/home/zyang/pg_workshop/4Sim_pgge/*.gfa
input_folder=/home/zyang/pg_workshop/4Sim_pgge
inputfa=/home/zyang/pg_workshop/4Sim.fa
output=/home/zyang/pg_workshop/4Sim_pgge
beehave=/home/zyang/pg_workshop/beehave.R






for x in $inputGFA

do
pgge -g $x -f $inputfa -o $output -r $beehave -b $output/pgge_4Sim_peanut_bed -l 100000 -s 5000 -t 16

done
