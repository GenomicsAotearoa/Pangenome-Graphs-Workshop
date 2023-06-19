# 6. Variant calling

## From graphs to VCFs


### `vg deconstruct` graph to get the variations in vcf 
<!-- 
!!! terminal "code"

    ```bash
    mkdir vg_deconstruct
    #copy gfa to vg_deconsturct and rename the gfa files with its PGGB settings
    cp /home/zyang/pg_workshop/4Sim_1K96/4Sim.fa.97e7156.417fcdf.7659dc8.smooth.final.gfa /home/zyang/pg_workshop/vg_deconstruct/4Sim_1K96.gfa
    cp /home/zyang/pg_workshop/4Sim_1K96_K79/4Sim.fa.f958389.417fcdf.7659dc8.smooth.final.gfa /home/zyang/pg_workshop/vg_deconstruct/4Sim_1K96_K79.gfa
    
    cd /home/zyang/pg_workshop/vg_deconstruct
    
    module purge
    module load vg/1.46.0
    module load BCFtools/1.15.1-GCC-11.3.0
    
    vg deconstruct -p NC_017518  -a -e 4Sim_1K96.gfa >4Sim_1K96_aep1.vcf
    bcftools stats 4Sim_1K96_aep1.vcf >4Sim_1K96_aep1.vcf_stats
    ```
-->

Set up directory for VG and GFA files for each pangenome graph.

!!! terminal "code"

    ```bash
    # Return to working directory
    cd ~/pg_workshop

    # Create VG directory
    mkdir -p vg_deconstruct
    ```

Load the necessary modules for an example run.

!!! terminal "code"

    ```bash
    module purge
    module load vg/1.46.0
    module load BCFtools/1.15.1-GCC-11.3.0
    ```

An example run to obtain VCF files from GFA.

!!! terminal "code"

    ```bash
    vg deconstruct -p NC_017518 -a -e 4Sim_1K96/4Sim.fa.97e7156.417fcdf.7659dc8.smooth.final.gfa > 4Sim_1K96_aep1.vcf
    bcftools stats 4Sim_1K96_aep1.vcf > 4Sim_1K96_aep1.vcf_stats
    ```
