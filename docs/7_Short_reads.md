# 7. Short reads
NGS data analysis used graph as a reference 

## vg mapping preliminaries
Although vg contains a number of tools for working with pangenome graphs, it is best-known for read mapping. This is ultimately what many of its users are interested in `vg` for. In fact, vg contains three mature short read mapping tools:

!!! info ""

    - `vg map`: the original, highly accurate mapping algorithm
    - `vg giraffe`: the much faster and still accurate haplotype-based mapping algorithm
    - `vg mpmap`: the splice-aware RNA-seq mapping algorithm

more details of vg can be found https://github.com/vgteam/vg

we use vg map in this workshop 

### Learning objectives

!!! quote ""

    - Map NGS data to graph using vg map
    - Variant calling for NGS data against genome graph 


## Build index for graph

!!! terminal "code"

    ```bash
    
    mkdir graph_NGS
    
    #copy graph to the graph reference (.gfa file) to work direvtory graph_NGS 
    cp /home/$your_home_dir/pg_workshop/5NM*.gfa ./home/$your_home_dir/pg_workshop/graph_NGS/5NM.gfa

    cd /home/$your_home_dir/pg_workshop/graph_NGS
    ```

Load the necessary modules for an example run.
!!! terminal "code"

    ```bash
    module purge
    module load vg/1.46.0
    ```


Build the index.

!!! terminal "code"

    ```bash
   
    mkdir -p ${temp_dir}
    
    # Convert graph into 256 bp chunks, saving as vg format
    vg mod -X 256 5NM.gfa > 5NM_256.vg

    #small graph is ok without prunning, complex graph will need to prune first before generating index
    #Build xg and gcsa index
    vg index -b ${temp_dir} -t 4 -x 5NM_256.xg -g 5NM_256.gcsa -k 16 5NM_256.vg
    ### you may have run out of temporary disk space at temp_dir

    ### pruning: use -M if pruning fails
    vg prune -u -m node-mapping.tmp -t 4 -k 24 5NM_256.vg > 5NM_256_chopped.vg
    
    vg index 5NM_256_chopped.vg -x 5NM_256_chopped.xg
    ### gcsa index, it takes .......
    
    vg index -b temp_dir -g 5NM_256_chopped.gcsa -x 5NM_256_chopped.xg -g 5NM_256_chopped.gcsa -k 16 5NM_256_chopped.vg
    ```

??? terminal "code"

    ```bash
    #run a slurm job for this, it takes ~10 mins based on the following setting 

    
    #!/bin/bash

    #SBATCH --account       nesi02659
    #SBATCH --job-name      build_index_for_5NMGraph
    #SBATCH --cpus-per-task 8
    #SBATCH --mem           16G
    #SBATCH --time          1:00:00
    #SBATCH --error         %x_%j.err
    #SBATCH --output        %x_%j.out

    # Modules
    module purge
    module load vg/1.46.0

    # Variables
    #cd /home/zyang/pg_workshop/graph_NGS
    data=5NM.gfa
    mkdir -p temp_dir

    # Convert graph into 256 bp chunks, saving as vg format
    vg mod -X 256 5NM.gfa > 5NM_256.vg

    #small graph is ok without prunning
    # Build xg and gcsa index
    #vg index -b ${temp_dir} -t $cpus-per-task -x 5NM_256.xg -g 5NM_256.gcsa -k 16 5NM_256.vg

    #complex graph will need to prune first before generating index
    ### pruning: use -M if pruning fails
    vg prune -u -m node-mapping.tmp -t 8 -k 24 5NM_256.vg > 5NM_256_chopped.vg

    vg index 5NM_256_chopped.vg -x 5NM_256_chopped.xg
    ### gcsa index
    vg index -b temp_dir -t 8 -x 5NM_256_chopped.xg -g 5NM_256_chopped.gcsa -k 16 5NM_256_chopped.vg

    ```



## Map NGS reads to graph 

Map reads back to graph reference
!!! terminal "code"

    ```bash
    # Modules
    module purge
    module load vg/1.46.0


     # Map reads
    vg map -t 8 -d 5NM_256_chopped -f NMI138_S5_R1_P.fastq.gz -f NMI138_S5_R1_P.fastq.gz -N NMI138  > NM138.vgmap_5NM.gam

    # Output mapping statistics
    vg stats -a NM138.vgmap_5NM.gam > /NM138.vgmap_4Sim_stats 
 
    ```

Map reads back to graph reference as a slurm job.
??? terminal "code"

    ```bash
    #!/bin/bash -e

    #SBATCH --account       nesi02659
    #SBATCH --job-name      vgmap_5e_5NM
    #SBATCH --cpus-per-task 24
    #SBATCH --mem           4G
    #SBATCH --time          01:00:00
    #SBATCH --error         %x_%j.err
    #SBATCH --output        %x_%j.out

    # Modules
    module purge
    module load vg/1.46.0

    # Variables
    wkdir=~/pg_test/graph_NGS
    index=${wkdir}/5NM_256_chopped.gcsa
    index_prefix=${index%%.gcsa

    # Map reads
    vg map -t $SLURM_CPUS_PER_TASK -d ${index_prefix} -f NMI138_S5_R1_P.fastq.gz -f NMI138_S5_R2_P.fastq.gz -N NMI138 > NM138.vgmap_5NM.gam

    # Output mapping statistics
    vg stats -a NM138.vgmap_5NM.gam > NM138.vgmap_5NM_stats
  
    ```
<!-- 3 min per read pair -->



## Genotying known variants 

!!! terminal "code"

    ```bash
    mkdir ~/pg_workshop/graph_NGS/vgmap_12e_sim4_allR10S3_typing
    ```

Generate snarls of graph.

!!! terminal "code"

    ```bash
    cd ~/pg_workshop/graph_NGS/refs

    module load vg/1.46.0

    vg snarls -t 2 4Sim_1K96_256.xg > 4Sim_1K96_256.xg.snarls

    cd ../
    ```

Perform genotyping.

!!! terminal "code"

    ```bash
    #!/bin/bash -e
    #SBATCH --account       nesi02659
    #SBATCH --job-name      5e_vgmap_genotying
    #SBATCH --cpus-per-task 24
    #SBATCH --mem           4G
    #SBATCH --time          01:00:00
    #SBATCH --error         %x_%j.err
    #SBATCH --output        %x_%j.out
    #SBATCH --array         0-5

    # Modules
    module purge 
    module load vg/1.46.0

    # Variables
    wkdir=~/pg_workshop/graph_NGS
    gam_dir=${wkdir}/graph_based_mapping
    out_dir=${wkdir}/vgmap_12e_sim4_allR10S3_typing

    mkdir -p ${out_dir}

    xg=${wkdir}/refs/4Sim_1K96_256.xg
    snarls=${wkdir}/refs/4Sim_1K96_256.xg.snarls

    # Array
    file_array=(${gam_dir}/*.gam)
    file=${file_array[$SLURM_ARRAY_TASK_ID]}
    prefix=$(basename ${file} .wgsim_er0.005.vgmap_4Sim.gam)

    # Calculate support reads
    vg pack -t $SLURM_CPUS_PER_TASK -x ${xg} -g ${file} -o ${out_dir}/${prefix}.vgmap_sim4_256_aln.pack

    # Call variants using the same coordinates and include reference calls for comparison
    vg call -t $SLURM_CPUS_PER_TASK -m 3,10 ${xg} -k ${out_dir}/${prefix}.vgmap_sim4_256_aln.pack -r $  {snarls} -a > ${out_dir}/${prefix}.vgmap_sim4_256_aln.pack_allR10S3.vcf
    ```

<!-- 3-5 min per input -->

<!-- 
```bash
#!/bin/bash

#SBATCH --account       nesi02659
#SBATCH --job-name      5e_vgmap_genotying
#SBATCH --cpus-per-task 8
#SBATCH --mem           4G
#SBATCH --time          24:00:00

module purge
module load vg/1.46.0

data_gam=/home/zyang/pg_workshop/graph_NGS/graph_based_mapping/*.wgsim_er0.005.vgmap_4Sim.gam
input=/home/zyang/pg_workshop/graph_NGS/graph_based_mapping
output=/home/zyang/pg_workshop/graph_NGS/vgmap_12e_sim4_allR10S3_typing
graph_xg=/home/zyang/pg_workshop/graph_NGS/refs/4Sim_1K96_256.xg
snarls_file=/home/zyang/pg_workshop/graph_NGS/refs/4Sim_1K96_256.xg.snarls



#compute snarls
#vg snarls $graph_xg >$snarls_file

for f in $data_gam
do

x=$(basename $f .wgsim_er0.005.vgmap_4Sim.gam)
echo ${x}


#Calculate the surpport reads ingoring mapping and base quality <5
#vg pack -t 48 -x $graph_xg -g $input/${x}.wgsim_er0.005.vgmap_4Sim.gam -Q 5 -o $output/${xvgmap_Sim4_256_aln.pack

#Calculate the surpport reads
vg pack -t 12 -x $graph_xg -g $input/${x}.wgsim_er0.005.vgmap_4Sim.gam -o $output/${xvgmap_sim4_256_aln.pack

#call variant using the same coordinates and including reference calls (for following compare)
vg call -t 12 -m 3,10 $graph_xg -k $output/${x}vgmap_sim4_256_aln.pack -r $snarls_file -a >$output/${x}vgmap_sim4_256_aln.pack_allR10S3.vcf

done  
```
-->

## Novel variant calling using graph reference

!!! terminal "code"

    ```bash
    mkdir /home/zyang/pg_workshop/graph_NGS/vgmap_5e_sim4_allR10S3_novelcalling 
    ```

!!! terminal "code"

    ```bash
    #!/bin/bash -e
    #SBATCH --account       nesi02659
    #SBATCH --job-name      5.call_novel_variant
    #SBATCH --cpus-per-task 24
    #SBATCH --mem           4G
    #SBATCH --time          01:00:00
    #SBATCH --error         %x_%j.err
    #SBATCH --output        %x_%j.out
    #SBATCH --array         0-5

    # Modules
    module purge 
    module load vg/1.46.0

    # Variables
    wkdir=~/pg_workshop/graph_NGS
    gam_dir=${wkdir}/graph_based_mapping
    out_dir=${wkdir}/vgmap_5e_sim4_allR10S3_novelcalling

    mkdir -p ${out_dir}

    vg=${wkdir}/refs/4Sim_1K96_256.vg
    xg=${wkdir}/refs/4Sim_1K96_256.xg

    # Array
    file_array=(${gam_dir}/*.gam)
    file=${file_array[$SLURM_ARRAY_TASK_ID]}
    prefix=$(basename ${file} .wgsim_er0.005.vgmap_4Sim.gam)

    # In order to also consider novel variants from the reads, use the augmented graph and gam 
    # (as created in the "Augmentation" example using vg augment -A).
    # Augment the graph with all variation from the GAM, saving to aug.vg
    ### Augment the graph with all variation form the GAM except 
    ### that implied by soft clips, saving to aug.vg.
    ### *aug-gam contains the same reads as aln.gam but mapped to aug.vg

    # Augment graph
    vg augment -t $SLURM_CPUS_PER_TASK ${vg} ${file} -A ${out_dir}/${prefix}.nofilt_aug.gam > ${out_dir}/${prefix}.nofilt_aug.vg

    # Index the augmented graph
    vg index -t $SLURM_CPUS_PER_TASK ${out_dir}/${prefix}.nofilt_aug.vg -x ${out_dir}/${prefix}.nofilt_aug.xg

    # Compute the all read support from the augmented GAM
    vg pack -t $SLURM_CPUS_PER_TASK -x ${out_dir}/${prefix}.nofilt_aug.xg -g ${out_dir}/${prefix}.nofilt_aug.gam -o ${out_dir}/${prefix}.nofilt_aug_allR.pack

    # Call variants
    vg call -t $SLURM_CPUS_PER_TASK -m 3,10 ${out_dir}/${prefix}.nofilt_aug.xg -k ${out_dir}/${prefix}.nofilt_aug_allR.pack > ${out_dir}/${prefix}.nofilt_aug_allR.pack.vcf
    ```

<!-- 16-18 min per gam -->

<!-- 
    ```bash
    #!/bin/bash
    
    #SBATCH --account       nesi02659
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
    vg pack -t 12 -x $output/${x}nofilt_aug.xg -g $output/${x}nofilt_aug.gam -o $output/${x}nofilt_aug_allR.pack
    
    
    #call variant
    vg call -t 12 -m 3,10 $output/${x}nofilt_aug.xg -k $output/${x}nofilt_aug_allR.pack >$output/${x}nofilt_aug_allR.pack.vcf
    
    #call variant snarl using the same coordinate
    #vg call -t 48 -m 3,10 $output/${x}nofilt_aug.xg -k $output/${x}nofilt_aug_allR.pack -a >$output/${x}nofilt_aug_allR.pack_snarls.vcf
    
    done
    ```
-->


