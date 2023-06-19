# PanGenome Graph Builder [(PGGB)](https://github.com/pangenome/pggb)
<p align="justify">
To investigate and analyse Neisseria Bacteria the PanGenome Graph Builder has been used. PanGenome Graph Builder (PGGB) is a computational tool used in genomics research to construct pan-genome graphs from large sets of genomic data. A pan-genome is the collection of all the genes and non-coding sequences present in a species or a group of related organisms. PGGB constructs a pan-genome graph, which is a data structure that represents the entire set of genes and genetic variations in a population. This graph can be used to study genetic diversity, gene function, and evolution. PGGB is designed to be scalable and efficient, making it suitable for large-scale genomic analyses. It is an open-source tool that can be used freely by researchers in the field of genomics.
</p>

### How does the pggb graph build work?
The pggb algorithm facilitates the construction of these graphs by progressively integrating genetic variants into a reference genome.
The key features and purpose of the PGGb include:
1.	Progressive Approach: The PGGb algorithm follows a step-by-step approach, iteratively adding genetic variants to the graph one at a time.
2.	Variant Integration: The tool efficiently integrates genetic variants, including single nucleotide variations (SNVs), insertions, deletions, and larger structural variations, into the graph representation.
3.	Optimal Placement: For each variant, the PGGb algorithm determines the most suitable position to incorporate it into the graph. This involves aligning the variant sequence with the existing graph structure while minimizing conflicts with the nodes and edges.
4.	Graph Expansion: Once the optimal placement is determined, the PGGb algorithm expands the graph by adding new nodes and edges to represent the variant sequence. The overall graph structure is modified to connect the variant nodes with the adjacent reference nodes.
5.	Large-Scale Graph Construction: The PGGb algorithm is designed to handle large-scale genomes and can efficiently construct genome graphs containing extensive genetic variations.

### All-to-all alignment
Generally, refers to the process of aligning all sequences in a given set against each other, rather than aligning them to a single reference sequence.
We begin with an alignment, with `wfmash`. This compares all sequences to each other and finds the best N mappings for each. It produces base-level alignments.

### Inducing the graph
Refers to the process of constructing the genome graph by progressively integrating genetic variants into a reference genome.
These base-level alignments are converted into a graph with `seqwish`. A filter is applied to remove short matches, which anchors the graph on confident longer exact matches.

### Normalizing the graph
refers to a process that aims to optimize the structure and representation of the genome graph by resolving redundant or overlapping elements. This step is typically performed after the initial construction of the graph.
The normalization process in PGGb involves several steps, which may vary depending on the specific implementation or version of the tool. Here are some common steps involved in normalizing the graph:
1.	Removal of Redundant Nodes: During the construction of the genome graph, it is possible that some nodes become redundant due to overlapping or repetitive sequences. Normalization involves identifying and removing these redundant nodes, streamlining the graph structure.
2.	Edge Optimization: In the graph, edges represent connections between nodes. During normalization, the edges are optimized to minimize redundancy and improve the efficiency of the graph. This can include merging or repositioning edges to create a more streamlined and accurate representation of the genome.
3.	Compact Representation: Normalization aims to reduce the overall size of the graph by compacting the representation. This can involve compressing repetitive regions or simplifying complex structures while preserving the essential information and variant representation.
4.	Graph Refinement: The normalization process also involves refining the graph structure by resolving inconsistencies, correcting errors, and improving the overall quality of the graph representation. This may include resolving conflicts between nodes and edges, addressing mismatches, and ensuring the graph accurately reflects the underlying genetic variations

To normalize the graph and harmonize the allele representation, we use `smoothxg` to apply a local MSA <!-- define: Multiple sequence alignment? --> across all parts of the graph.

### Downstream
These graphs offer a wide range of capabilities. Initially, we can generate several diagnostic visualizations derived from the graphs, providing a user-friendly way to comprehend the alignments at a broader level. Additionally, using the `PGGB` tool, we can generate variant calls by leveraging `vg deconstruct`. The graphs produced by `PGGB` serve as reference systems for aligning short reads through vg giraffe or long reads through `GraphAligner`. Furthermore, `odgi` allows us to utilize the graphs as reference systems to elucidate homology relationships across entire genomes.

### PanGenome Graph Evaluator [(PGGE)](https://github.com/pangenome/pgge)
<p align="justify">
The pangenome graph evaluation is a pipeline in the pangenome graphs, which measures the reconstruction accuracy of the graph. It helps find the best pangenome graph using input data and tasks.
</p>
