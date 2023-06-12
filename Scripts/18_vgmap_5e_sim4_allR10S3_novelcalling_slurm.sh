#!/usr/bin/bash

#SBATCH --account       ga03793
#SBATCH --job-name      5e_vgmap_novelvariant_calling 
#SBATCH --cpus-per-task 8
#SBATCH --mem           4G
#SBATCH --time          24:00:00

module purge
module load vg/1.46.0

data_gam=/home/zyang/pg_workshop/graph_NGS/graph_based_mapping/*.wgsim_er0.005.vgmap_4Sim.gam
input=/home/zyang/pg_workshop/graph_NGS/graph_based_mapping
output=/home/zyang/pg_workshop/graph_NGS/vgmap_5e_sim4_allR10S3_novelcalling
graph_vg=/home/zyang/pg_workshop/graph_NGS/refs/4Sim_1K96_256.vg
graph_xg=/home/zyang/pg_workshop/graph_NGS/refs/4Sim_1K96_256.xg


#compute snarls
#vg snarls $graph_xg >$output/${graph_xg}.snarls

for f in $data_gam
do 

x=$(basename $f .wgsim_er0.005.vgmap_4Sim.gam)
echo ${x}


#in order to also consider novel variants from the reads, use the augmented graph and gam (as created in the "Augmentation" example using vg augment -A)
#Augment augment the graph with all variation from the GAM, saving to aug.vg
### augment the graph with all variation from the GAM except 
### that implied by soft clips, saving to aug.vg
### *aug-gam contains the same reads as aln.gam but mapped to aug.vg

vg augment -t 12 $graph_vg $input/${x}.wgsim_er0.005.vgmap_4Sim.gam -A $output/${x}nofilt_aug.gam >$output/${x}nofilt_aug.vg

#index the augmented graph
vg index -t 12 $output/${x}nofilt_aug.vg -x $output/${x}nofilt_aug.xg

## Compute the all read support from the augmented gam 
vg pack -t 12 -x $output/${x}nofilt_aug.xg -g $output/${x}nofilt_aug.gam  -o $output/${x}nofilt_aug_allR.pack


#call variant
vg call -t 12 -m 3,10 $output/${x}nofilt_aug.xg -k $output/${x}nofilt_aug_allR.pack >$output/${x}nofilt_aug_allR.pack.vcf

#call variant snarl using the same coordinate
#vg call -t 48 -m 3,10 $output/${x}nofilt_aug.xg -k $output/${x}nofilt_aug_allR.pack -a >$output/${x}nofilt_aug_allR.pack_snarls.vcf

done 
