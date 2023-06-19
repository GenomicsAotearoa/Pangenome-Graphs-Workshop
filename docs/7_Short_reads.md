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

## Simulate NGS read data

We will use `wgsim` from SAMtools to simulate 2 $\times$ 150 bp NGS data (based on 5 genomes) with an error rate of 0.005.

!!! terminal "code"

    ```bash
    mkdir ~/pg_workshop/graph_NGS
    cd graph_NGS
    mkdir simu_NGS_data
    ```
The script for simulation NGS data

!!! terminal "code"

    === "Loop"

        ```bash
        #!/bin/bash -e
        #SBATCH --account       nesi02659
        #SBATCH --job-name      5.simulate_NGS
        #SBATCH --cpus-per-task 2
        #SBATCH --mem           4G
        #SBATCH --time          30:00
        #SBATCH --error         %x_%j.err
        #SBATCH --output        %x_%j.out

        # Modules
        module purge
        module load SAMtools/1.16.1-GCC-11.3.0

        # Variables
        wkdir=~/pg_workshop
        input_folder=${wkdir}/dataset_for_pg_workshop/12_genomes_for_NGS_simulation
        output_folder=${wkdir}/graph_NGS/simu_NGS_data

        mkdir -p ${output_folder}

        # Run
        for input in ${input_folder}/*_6k.fa; do
          # Set prefix
          prefix=$(basename ${input} .fa)
          # Simulate reads
          wgsim -N 1000000 -1 150 -2 150  -e 0.005 -r 0 -R 0 -X 0 \
            ${input} \
            ${output_folder}/${prefix}.wgsim_er0.005.R1.fq \
            ${output_folder}/${prefix}.wgsim_er0.005.R2.fq
          # Compress reads
          gzip ${output_folder}/${prefix}.wgsim_er0.005.R*.fq
        done
        ```

    === "Array"

        ```bash
        #!/bin/bash -e
        #SBATCH --account       nesi02659
        #SBATCH --job-name      5.simulate_NGS
        #SBATCH --cpus-per-task 2
        #SBATCH --mem           4G
        #SBATCH --time          05:00
        #SBATCH --error         %x_%j_%a.err
        #SBATCH --output        %x_%j_%a.out
        #SBATCH --array         0-5

        # Modules
        module purge
        module load SAMtools/1.16.1-GCC-11.3.0

        # Variables
        wkdir=~/pg_workshop
        input_folder=${wkdir}/dataset_for_pg_workshop/12_genomes_for_NGS_simulation
        output_folder=${wkdir}/graph_NGS/simu_NGS_data

        mkdir -p ${output_folder}

        # Array
        input_array=(${input_folder}/*_6k.fa)
        input=${input_array[$SLURM_ARRAY_TASK_ID]}
        prefix=$(basename ${input} .fa)

        # Simulate reads
        wgsim -N 1000000 -1 150 -2 150  -e 0.005 -r 0 -R 0 -X 0 \
          ${input} \
          ${output_folder}/${prefix}.wgsim_er0.005.R1.fq \
          ${output_folder}/${prefix}.wgsim_er0.005.R2.fq

        # Compress reads
        gzip ${output_folder}/${prefix}.wgsim_er0.005.R*.fq
        ```

## Build index for graph

Copy a reference (4Sim_1K96.gfa) into the NGS read directory.

!!! terminal "code"

    ```bash
    mkdir refs
    
    # copy graph to the refs work direvtory 
    cp /home/zyang/pg_workshop/vg_deconstruct/4Sim_1K96.gfa /home/zyang/pg_workshop/graph_NGS/refs
    
    #make tem_dir
    mkdir /home/zyang/pg_workshop/graph_NGS/refs/temp_dir
    ```

Build the index.

!!! terminal "code"

    ```bash
    #!/bin/bash
    
    #SBATCH --account       nesi02659
    #SBATCH --job-name      build_index_for_4SimGraph
    #SBATCH --cpus-per-task 8
    #SBATCH --mem           4G
    #SBATCH --time          1:00:00
    #SBATCH --error         %x_%j.err
    #SBATCH --output        %x_%j.out
    
    # Modules
    module purge
    module load vg/1.46.0
    
    # Variables
    cd ~/pg_workshop/graph_NGS/refs
    data=4Sim_1K96.gfa
    temp_dir=temp_dir
    prefix=$(basename ${data} .gfa)

    mkdir -p ${temp_dir}

    # Convert graph into 256 bp chunks, saving as vg format
    vg mod -X 256 ${prefix}.gfa > ${prefix}_256.vg

    # Build xg and gcsa index
    vg index -b ${temp_dir} -t $SLURM_CPUS_PER_TASK -x ${prefix}_256.xg -g ${prefix}_256.gcsa -k 16 ${prefix}_256.vg
    
    #small graph is ok without prunning, complex graph will need to prune first before generating index
    
    ### pruning: use -M if pruning fails
    #vg prune -u -m node-mapping.tmp -t 48 -k 24 ${x}_256.vg > ${x}_256_chopped.vg
    
    #vg index ${x}_256_chopped.vg -x ${x}_256_chopped.xg
    ### gcsa index
    #vg index -b $tem_dir -t 48  -g ${x}_256_chopped.gcsa  ${x}_256_chopped.vg
    ```

## Map NGS reads to graph 

!!! terminal "code"

    ```bash
    mkdir /home/zyang/pg_workshop/graph_NGS/graph_based_mapping
    ```

Map reads back to reference pangenome graph as an array job.

!!! terminal "code"

    ```bash
    #!/bin/bash -e
    #SBATCH --account       nesi02659
    #SBATCH --job-name      vgmap_5e_4Sim
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
    read_folder=${wkdir}/simu_NGS_data
    output_folder=${wkdir}/graph_based_mapping
    index=${wkdir}/refs/4Sim_1K96_256.gcsa
    index_prefix=${index%%.gcsa}

    mkdir -p ${output_folder}

    # Array
    file_array=(${read_folder}/*R1.fq.gz)
    file=${file_array[$SLURM_ARRAY_TASK_ID]}
    read_prefix=$(basename ${file} .R1.fq.gz)
    read2=$(echo ${file} | sed -e 's/R1.fq.gz/R2.fq.gz/')

    # Map reads
    vg map -t $SLURM_CPUS_PER_TASK -d ${index_prefix} -f ${file} -f ${read2} -N ${read_prefix} > ${output_folder}/${read_prefix}.vgmap_4Sim.gam

    # Output mapping statistics
    vg stats -a ${output_folder}/${read_prefix}.vgmap_4Sim.gam > ${output_folder}/${read_prefix}.vgmap_4Sim_stats
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


