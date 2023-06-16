#!/bin/bash


WD=/nesi/nobackup/nesi02659/pg_workshop #Working Directory

data=${WD}/4Sim.fa

 


pgge -g ${WD}/output/*.gfa -f $data -o pgge -r ${WD}/beehave.R -b pgge/pgge_4Sim_peanut_bed -l 100000 -s 5000 -t 16
