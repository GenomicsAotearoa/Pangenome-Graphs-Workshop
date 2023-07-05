# 5. Evaluating output
When we ran PGGB, the 'odgi stats -S' option was utilized to generate statistics for both the seqwish and smoothxg graphs and the 'multiqc -m' option was employed to generate a MultiQC report, providing comprehensive statistics and visualizations of the graphs. All pertinent results can be located in the MultiQC report, which is saved in HTML format.


## Pangenome graph visualization using ODGI 

### ODGI Compressed 1D visualization
!!! info ""
  
![ODGI Compressed 1D visualization](theme_figures/ODGI-Compressed-1D-5NM.png)

This image shows a 1D rendering of the built pangenome graph. The graph nodes are arranged from left to right, forming the pangenome sequence. Summarization of path coverage across all paths. Dark blue means highest coverage. Dark red means lowest coverage. The path names are placed on the left. The black lines under the paths are the links, which represent the graph topology.

### ODGI 1D visualization
!!! info ""
 
![ODGI 1D visualization](theme_figures/ODGI-1D-5NM.png)

This image shows a 1D rendering of the built pangenome graph. The graph nodes are arranged from left to right, forming the pangenome sequence. The colored bars represent the paths versus the pangenome sequence in a binary matrix. The path names are placed on the left. The black lines under the paths are the links, which represent the graph topology.


### ODGI 1D visualization by path position
!!! info ""


![ODGI 1D visualization by path position](theme_figures/ODGI-Path-Position-1D-5NM.png)

This shows a 1D rendering of the built pangenome graph where the paths are colored according to their nucleotide position. Light grey means a low path position, black is the highest path position.

### ODGI 1D visualization by path orientation
!!! info ""

![ODGI 1D visualization by path orientation](theme_figures/ODGI-Path-Orientation-1D-5NM.png)
This image shows a 1D rendering of the built pangenome graph where the paths are colored by orientation. Forward is black, reverse is red.

 
??? info "What makes the last path different compared to the other paths?"

    The orientation of the last path is almost exactly the reverse of the second to last one, right? Do you think it's possible that the last path of the genome was submitted as its reverse complement? 



## Circlator

**As bacterial genomes are circular, if we can fix the start of the input genomes for pangenome graph constrcution, this may help to exclude the unneccessary complexity in the grph **
??? info ""

    let's fix the start for all genome using circlator, submit a slurm job. It takes less than one minute for each sample. 
    ```bash
    #!/usr/bin/bash

    #SBATCH --account       ga03793
    #SBATCH --job-name      restart_fna
    #SBATCH --cpus-per-task 8
    #SBATCH --mem           4G
    #SBATCH --time          1:00:00

    module load Circlator/1.5.5-gimkl-2022a-Python-3.10.5

    cd /home/zyang/pg_test
    data=/home/zyang/pg_test/*.fna

   for f in $data
   do

   x=$(basename $f .fna)
   echo ${x}

   circlator fixstart  ${x}.fna  ${x}.restart

  done
  ```

??? info "rebuild the "




    

### 1D visualization by node depth
!!! info ""

![ODGI 1D visualization by node depth](theme_figures/ODGI-Node-Depth-1D-5NM.png)
This shows a 1D rendering of the built pangenome graph where the paths are colored according to path depth. Using the Spectra color palette with 4 levels of path depths, white indicates no depth, while grey, red, and yellow indicate depth 1, 2, and greater than or equal to 3, respectively.


### ODGI 1D visualization by uncalled bases
!!! info ""

![ODGI 1D visualization by uncalled bases](theme_figures/ODGI-Uncalled-1D-5NM.png)
This shows a 1D rendering of the built pangenome graph where the paths are colored according to the coverage of uncalled bases. The lighter the green, the higher the 'N' content of a node is.


### ODGI 2D drawing
!!! info ""

![ODGI 2D visualization](theme_figures/ODGI-2D-5NM-small.png)



### Check the statistics statistics for both the seqwish and smoothxg graphs
!!! info ""

#### 5NM 2k94
| Sample Name                         | Length    | Nodes  | Edges  | Components | A   |C    |T    |G    |N   |
|:-----                               |----------:|-------:|-------:|-----------:|----:|----:|----:|----:|----:|
|seqwish	|2280344	|22216	|29823	|4	|1	|551639	|578590	|557450	|592665	|0|
|smooth	|2261163	|29965	|40179	|4	|1	|548754	|574693	|551650	|586066	|0|

![1k96 ODGI 1D visualization by path orientation](theme_figures/4Sim.fa.97e7156.417fcdf.7659dc8.smooth.final.og.viz_inv_multiqc.png)


#### 1k96,-K79


| Sample Name                         | Length    | Nodes  | Edges  |Paths       |Components | A   |C    |T    |G    |N   |
|:-----                               |----------:|-------:|-------:|------------|-----------:|----:|----:|----:|----:|----:|
|seqwish	|3165112	|123218	|166072	|5	|1	|786941	|789388	|784867	|803816	|100|
|smooth	  |2926040	|246217	|332089	|5	|1	|735148	|742974	|726837	|720981 |100|





## Circlator

**NB: this material below was moved from Section 2. Motivate this stuff by inspecting the 1D graph of inversions, and
note that one of the samples is completely inverted relative to the others.**

!!! question "Exercises"

    - Can we cat the graph and start pangenome graph construction now? 
    - What potiential issues could there be? 
    - We need to check whether all the genomes are with the same start. If not, it will cause unwanted complexsity for the pangenome graph. 
    
```bash
head -10 *.fna >>head10_check
less -S head10_check
```
let's fix the start for all genome using circlator, submit a slurm job. It takes less than one minute for each sample. 
```bash
#!/usr/bin/bash

#SBATCH --account       ga03793
#SBATCH --job-name      restart_fna
#SBATCH --cpus-per-task 8
#SBATCH --mem           4G
#SBATCH --time          1:00:00

module load Circlator/1.5.5-gimkl-2022a-Python-3.10.5

cd /home/zyang/pg_test
data=/home/zyang/pg_test/*.fna

for f in $data
do

x=$(basename $f .fna)
echo ${x}

circlator fixstart  ${x}.fna  ${x}.restart

done
```

