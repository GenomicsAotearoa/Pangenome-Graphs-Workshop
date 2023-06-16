#!/bin/bash


WD=/nesi/nobackup/nesi02659/pg_workshop #Working Directory

data=${WD}/ASM19152v1_pgsim.fa

 


singularity exec ${container} pggb -i $data -s 1000 -p 95 -n 6 -k 79 -t 2 -S -m -o output -V 'NC_017518.1:#' 
