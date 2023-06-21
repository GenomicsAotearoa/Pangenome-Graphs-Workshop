# 1. Overview of pangenome graphs 
<p align="justify">
A pangenome is defined as the comprehensive collection of whole-genome sequences from multiple individuals within a clade, a population or a species. This collective genomic dataset can be further divided into two distinct components: the core genome, which includes genes present in all individuals at the time of analysis, and the accessory genome, consisting of genes found only in a subset of individuals. 
</p>

![image](https://github.com/GenomicsAotearoa/Pangenome-Graphs-Workshop/blob/main/docs/theme_figures/Fig.1_pangenome.png)

<center><small>Overview of the pangenome graph workflow</small></center>


<p align="justify">
Pangenome graphs are pangenomes stored in graph models that can capture the entire genetic variation among genomes in a population or of a set of related organisms (Figure 1B). There are three components of a variation graph: Nodes, edges and paths.
</p>

#### **Nodes**
- DNA segments, which can be any length 


#### **Edges** 
- Describe the possible ways of walking through the nodes
- Connect pairs of node strands
- Can represent inversions 


#### **Paths** 
- Paths are routes through the nodes of the graph
- Genomes
- Haplotypes
- Alleles/variants 

<p align="justify">
This pipeline for pangenome graphs comprises three key stages: graph construction using `PGGB`, graph manipulation via `ODGI`, and variant calling for Next-Generation Sequencing (NGS) data utilizing the VG toolkit, as shown in Figure 1C.

The PGGB pipeline, which operates without a reference method, builds pangenome graphs using an all-to-all whole genome alignment approach with `wfmash`. Subsequent graph induction is accomplished through `seqwish`, followed by progressive normalization implemented with `smoothxg` and `gfaffix`.

ODGI is employed for various graph manipulation tasks, including visualization and the extraction of distances between paths within the graph. This feature enables further phylogenetic analysis.

By using the pangenome graph created with PGGB, it is possible to concurrently identify a variety of genetic variations. These include structural variations (SVs), rearrangements, and smaller variants such as single nucleotide polymorphisms (SNPs) and insertions/deletions, which can be identified through the `vg deconstruction` process.

Lastly, the VG toolkit is harnessed for NGS data analysis against the graph, which includes tasks like read mapping and variant calling.

</p>

![image](./theme_figures/Fig.1_overview%20of%20pangenome%20graph%20pipeline.png)

<center><small>Overview of the pangenome graph workflow</small></center>
