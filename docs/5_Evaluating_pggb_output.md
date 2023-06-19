# 5. Evaluating output
### Evaluate Pangenome Graphs for 4Sim Genomes Constructed with Different Settings
- We have employed three distinct settings to construct the pangenome graph of the 4Sim genomes. Which setting yielded the most optimal result? How can we determine this? 

- Download the multiqc.html file, check the Detailed ODGI stats table.
#### 1k96
| Sample Name                         | Length    | Nodes  | Edges  | Components | A   |C    |T    |G    |N   |
|:-----                               |----------:|-------:|-------:|-----------:|----:|----:|----:|----:|----:|
|seqwish	|2280344	|22216	|29823	|4	|1	|551639	|578590	|557450	|592665	|0|
|smooth	|2261163	|29965	|40179	|4	|1	|548754	|574693	|551650	|586066	|0|

![1k96 ODGI 1D visualization by path orientation](theme_figures/4Sim.fa.97e7156.417fcdf.7659dc8.smooth.final.og.viz_inv_multiqc.png)


#### 1k96,-K79


| Sample Name                         | Length    | Nodes  | Edges  | Components | A   |C    |T    |G    |N   |
|:-----                               |----------:|-------:|-------:|-----------:|----:|----:|----:|----:|----:|
|seqwish	|2279905	|22209	|29812	|4	|1	|551588	|578459	|557291	|592567	|0|
|smooth	|2261401	|29976	|40199	|4	|1	|553049	|580469	|547460	|580423	|0|


![1k96K79 ODGI 1D visualization by path orientation](theme_figures/4Sim.fa.f958389.417fcdf.7659dc8.smooth.final.og.viz_inv_multiqc.png)

#### 10k96
| Sample Name                         | Length    | Nodes  | Edges  | Components | A   |C    |T    |G    |N   |
|:-----                               |----------:|-------:|-------:|-----------:|----:|----:|----:|----:|----:|
|seqwish	|2340700	|22166	|29759	|4	|1	|566559	|594741	|571836	|607564	|0|
|smooth	|2319601	|29888	|40070	|4	|1	|566755	|599287	|562078	|591481	|0|

![10k96 ODGI 1D visualization by path orientation](theme_figures/4Sim.fa.e7f7fe6.417fcdf.7659dc8.smooth.final.og.viz_inv_multiqc.png)


