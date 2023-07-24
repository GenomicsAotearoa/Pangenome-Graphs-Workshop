# Short reads
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
    #SAMN13450731 is the NCBI record for NMI138 
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

Generate snarls of graph.

!!! terminal "code"

    ```bash
    module purge
    module load vg/1.46.0

    vg snarls -t 2 5NM_256.xg > 5NM_256.xg.snarls
    
    ```

Perform genotyping.


!!! terminal "code"

    ```bash
    # Calculate support reads
    vg pack -t 8 -x 5NM_256.xg -g NM138.vgmap_5NM.gam -o NM138_vgmap_5NM_256.pack

    # Call variants using the same coordinates and include reference calls for comparison
    vg call -t 8 -m 3,10 5NM_256.xg -k NM138_vgmap_5NM_256.pack -r 5NM_256.xg.snarls -a > NM138.vgmap_5NM_256.pack_allR10S3.vcf
    
    ```



??? terminal "code"

    ```bash
    #!/bin/bash -e
    #SBATCH --account       nesi02659
    #SBATCH --job-name      5NM_vgmap_genotying
    #SBATCH --cpus-per-task 24
    #SBATCH --mem           4G
    #SBATCH --time          01:00:00
    #SBATCH --error         %x_%j.err
    #SBATCH --output        %x_%j.out

    # Modules
    module purge 
    module load vg/1.46.0

    vg index 5NM_256.vg -x 5NM_256.xg

    vg snarls -t 2 5NM_256.xg > 5NM_256.xg.snarls

    # Calculate support reads
    vg pack -t 8 -x 5NM_256.xg -g NM138.vgmap_5NM.gam -o NM138_vgmap_5NM_256.pack

    # Call variants using the same coordinates and include reference calls for comparison
    vg call -t 8 -m 3,10 5NM_256.xg -k NM138_vgmap_5NM_256.pack -r 5NM_256.xg.snarls -a > NM138.vgmap_5NM_256.pack_allR10S3.vcf
  
    ```


## Novel variant calling using graph reference

!!! terminal "code"

    ```bash
     # In order to also consider novel variants from the reads, use the augmented graph and gam 
    # (as created in the "Augmentation" example using vg augment -A).
    # Augment the graph with all variation from the GAM, saving to aug.vg
    ### Augment the graph with all variation form the GAM except 
    ### that implied by soft clips, saving to aug.vg.
    ### *aug-gam contains the same reads as aln.gam but mapped to aug.vg

    # Augment graph
    vg augment -t 8  5NM_256_chopped.vg NM138.vgmap_5NM.gam -A NM138.nofilt_aug.gam > NM138.nofilt_aug.vg

    # Index the augmented graph
    vg index -t 8  NM138.nofilt_aug.vg  -x NM138.nofilt_aug.xg

    # Compute the all read support from the augmented GAM
    vg pack -t 8 -x NM138.nofilt_aug.xg -g NM138.nofilt_aug.gam -o NM138.nofilt_aug_allR.pack
    
    # Call variants.
    #we need to trouble shooting about this, why the vcf file is empty 
    vg call -t 8 -m 3,10 NM138.nofilt_aug.xg -k NM138.nofilt_aug_allR.pack > NM138.nofilt_aug_allR.pack.vcf
    
    ```

??? terminal-2 "Slurm script"

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



