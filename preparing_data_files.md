# Preparing required data files with simuG tool
[simuG](https://github.com/yjx1217/simuG) is a light-weighted tool which can be used to simulate sequences with pre-defined number of random genomic variants against a reference sequence (Please refer the [Git Page](https://github.com/yjx1217/simuG) for more details). The following example demostrates generating three new sequence using known SNP and INDEL counts.

### 1. Setting up the tool 
Please refer the installation instruction in the [Git Page](https://github.com/yjx1217/simuG).
```bash
git clone https://github.com/yjx1217/simuG.git
cd simuG
perl simuG.pl -h
perl vcf2model.pl -h
```

### 2. Data source
You can download the reference sequence (Neisseria gonorrhoeae FA 1090 genome assembly) from the site [Genome assembly ASM19152v1](https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_000191525.1/) and get the file _ncbi_dataset/data/GCF_000191525.1/GCF_000191525.1_ASM19152v1_genomic.fna_ from he zip file GCF_000191525.1.zip. Rename it as "ref.fa" to refere easily in next steps. 

In Unix environment you can use `curl`. 

```bash
$ curl -OJX GET "https://api.ncbi.nlm.nih.gov/datasets/v2alpha/genome/accession/GCF_000191525.1/download?include_annotation_type=GENOME_FASTA,GENOME_GFF,RNA_FASTA,CDS_FASTA,PROT_FASTA,SEQUENCE_REPORT&filename=GCF_000191525.1.zip" -H "Accept: application/zip"
  % Total    % Received % Xferd  Average Speed   Time    Time     Time  Current
                                 Dload  Upload   Total   Spent    Left  Speed
100 1807k    0 1807k    0     0   294k      0 --:--:--  0:00:06 --:--:--  421k

$ ls
GCF_000191525.1.zip

$ unzip GCF_000191525.1.zip 
Archive:  GCF_000191525.1.zip
  inflating: README.md               
  inflating: ncbi_dataset/data/assembly_data_report.jsonl  
  inflating: ncbi_dataset/data/GCF_000191525.1/GCF_000191525.1_ASM19152v1_genomic.fna  
  inflating: ncbi_dataset/data/GCF_000191525.1/genomic.gff  
  inflating: ncbi_dataset/data/GCF_000191525.1/cds_from_genomic.fna  
  inflating: ncbi_dataset/data/GCF_000191525.1/protein.faa  
  inflating: ncbi_dataset/data/GCF_000191525.1/sequence_report.jsonl  
  inflating: ncbi_dataset/data/dataset_catalog.json  

$ cp ncbi_dataset/data/GCF_000191525.1/GCF_000191525.1_ASM19152v1_genomic.fna ./

$ head GCF_000191525.1_ASM19152v1_genomic.fna 
>NC_017518.1 Neisseria meningitidis NZ-05/33, complete sequence
TTCGGCTTAAACCTTATCCATATCCAAACGCATAACCGTAACCCATTCACCGTTATGGAAATGTCGCCCGACAACCGCCC
AGCCGAATGATTCATAAAATATTTGCACATCAGGCGTATAAAGATACAAGAACTTTATCCCCAGCGAACGCGCTGCGCCT
ATGCAGTGGGCGACCAGCCTCCTGCCAATGCCTTTTCCGCGATATTCAGGTAAAACAAAGACATCACCCAACCAATATTC
ATACTGTGGAAAACTTTCCATATCATGCCGCTTGACCGTAGCCGAACCCAACAGGGTTCCGGAATCATCCACAGCCGCAA
AAGCCAGCGGCAGTTCGTCATCCTTCAAACACCTGCCGTAATAGGCATGAATCTTATCCACAGAAGACCACGGTTCAAAT
CCGTGCCACTCCTCAAACAACGCCTGAACCAACCTGCCGATATGCCCGGCTTTCAGCCGTGTAATGAAAACAGTATTGTC
CACAAAGAGGGAATTCATCGGTCAATTCCCCGACGCTTTCGTTCCCCCTGCGCCGTAAACCGCATTCCAAGCATAGTCCA
AACGCACTCCGATTTGCCTCAGCTCTTCAGCCTGCCGGGCTTTTTGCGCCATTGCTGCAGGAATTTCCGCTTCCAAACGG
GCGATGTCTGCCTGAGCCGTCTGCAAACGCCGGCGCGCATCTTCCAAATCCGACTGCATCCCGATGATTTTTCCGTCCAG
```

### 3. Simulating sequences with known genetics variants counts using simuG
#### i. SNP=5000 
```
$ perl simuG.pl -refseq GCF_000191525.1_ASM19152v1_genomic.fna -snp_count 5000 -prefix Simulation_SNP_5000 -seed 112345678900

$ ls -1sh Simulation_SNP_5000*
512K Simulation_SNP_5000.refseq2simseq.map.txt
768K Simulation_SNP_5000.refseq2simseq.SNP.vcf
2.3M Simulation_SNP_5000.simseq.genome.fa

$ less Simulation_SNP_5000.simseq.genome.fa 
>NC_017518.1
TTCGGCTTAAACCTTATCCATATCCAAACGCATAACCGTAACCCATTCACCGTTATGGAAATGTCGCCCGACAACCGCCCAGCCGAATGATTCATAAAATATTTGCACATCAGGCGTATAAAGATACAAGAACTTTATCCCCAGCGAACGCGCTGCGCCTATGCAGTGGGCGACCAGCCTCCTGCCAATGCCTTTTCCGCGATATTCAGGTAAAACAAAGACATCAC
```
After the simulation the Chromosome number is same. Let's rename it be more clear in the next steps. We can use `sed -i` for this. 
```
$ sed -i '/>/ s/\(.*\)/\1_SNP_5000/' Simulation_SNP_5000.simseq.genome.fa

$ less Simulation_SNP_5000.simseq.genome.fa
>NC_017518.1_SNP_5000
TTCGGCTTAAACCTTATCCATATCCAAACGCATAACCGTAACCCATTCACCGTTATGGAAATGTCGCCCGACAACCGCCCAGCCGAATGATTCATAAAATATTTGCACATCAGGCGTATAAAGATACAAGAACTTTATCCCCAGCGAACGCGCTGCGCCTATGCAGTGGGCGACCAGCCTCCTGCCAATGCCTTTTCCGCGATATTCAGGTAAAACAAAGACATCAC
```

#### ii. INDEL=5000 (IN:DEL ratio 1:1)
```bash
$ perl simuG.pl -refseq GCF_000191525.1_ASM19152v1_genomic.fna -indel_count 5000 -prefix Simulation_INDEL_5000 -seed 212345678900

$ ls -1tr Simulation_INDEL_5000*
Simulation_INDEL_5000.refseq2simseq.INDEL.vcf
Simulation_INDEL_5000.refseq2simseq.map.txt
Simulation_INDEL_5000.simseq.genome.fa

$ sed -i '/>/ s/\(.*\)/\1_INDEL_5000/' Simulation_INDEL_5000.simseq.genome.fa

$ less Simulation_INDEL_5000.simseq.genome.fa
>NC_017518.1_INDEL_5000
TTCGGCTTAAACCTTATCCATATCCAAACGCATAACCGTAACCCATTCACCGTTATGGAAATGTCGCCCGACAACCGCCCAGCCGAATGATTCATAAAATATTTGCACATCAGGCGTATAAAGATACAAGAACTTTATCCCCAGCGAACGCGCTGCGCCTATGCAGTGGGCGACCAGCCTCCTGCCAATGCCTTTTCCGCGATATTCAGGTAAAACAAAGACATCACC
```
_note: Remeber to use different seeds for each simulation_

#### iii. SNP=4000 and INDEL=4000 (IN:DEL ratio 1:4)
```bash
$ perl simuG.pl -refseq GCF_000191525.1_ASM19152v1_genomic.fna -snp_count 4000 -indel_count 4000 -ins_del_ratio 0.25 -prefix Simulation_SNP_4000_INDEL_4000 -seed 312345678900

$ ls -1tr Simulation_SNP_4000_INDEL_4000*
Simulation_SNP_4000_INDEL_4000.refseq2simseq.SNP.vcf
Simulation_SNP_4000_INDEL_4000.refseq2simseq.INDEL.vcf
Simulation_SNP_4000_INDEL_4000.refseq2simseq.map.txt
Simulation_SNP_4000_INDEL_4000.simseq.genome.fa

$ sed -i '/>/ s/\(.*\)/\1_SNP_4000_INDEL_4000/' Simulation_SNP_4000_INDEL_4000.simseq.genome.fa

$ less Simulation_SNP_4000_INDEL_4000.simseq.genome.fa
>NC_017518.1_SNP_4000_INDEL_4000
TTCGGCTTAAACCTTATCCATATCCAAACGCATAACCGTAACCCATTCACCGTTATGGAAATGTCGCCCGACAACCGCCCAGCCGAATGATTCATAAAATATTTGCACATCAGGCGTATAAAGATACAAGAACTTTATCCCCAGCGAACGCGCTGCGCCTATGCAGTGGGCGACCAGCCTCCTGCCAATGCCTTGTCCGCGATATTCAGGTAAAACAAAGACATCACCCAAC
```

#### iv. SNP=4000, INDEL=4000 (IN:DEL ratio 1:4) and 4 inversions
To have a more complex sample we'll apply 4 inversions to the above simulated sequence. 
```bash
$ perl simuG.pl -refseq Simulation_SNP_4000_INDEL_4000.simseq.genome.fa -inversion_count 4 -inversion_min_size 50000 -prefix Simulation_SNP_4000_INDEL_4000_INV_4 -seed 412345678900

$ ls -1tr Simulation_SNP_4000_INDEL_4000_INV_4*
Simulation_SNP_4000_INDEL_4000_INV_4.refseq2simseq.inversion.vcf
Simulation_SNP_4000_INDEL_4000_INV_4.refseq2simseq.map.txt
Simulation_SNP_4000_INDEL_4000_INV_4.simseq.genome.fa

$ sed -i '/>/ s/\(.*\)/\1_INV_4/' Simulation_SNP_4000_INDEL_4000_INV_4.simseq.genome.fa

$ less Simulation_SNP_4000_INDEL_4000_INV_4.simseq.genome.fa
>NC_017518.1_SNP_4000_INDEL_4000_INV_4
TTCGGCTTAAACCTTATCCATATCCAAACGCATAACCGTAACCCATTCACCGTTATGGAAATGTCGCCCGACAACCGCCCAGCCGAATGATTCATAAAATATTTGCACATCAGGCGTATAAAGATACAAGAACTTTATCCCCAGCGAACGCGCTGCGCCTATGCAGTGGGCGACCAGCCTCCTGCCAATGCCTTGTCCGCGATATTCAGGTAAAACAAAGACATCAC
```

#### iv. SNP=4000, INDEL=4000 (IN:DEL ratio 1:4) and 4 copy number variations 
Let's apply 4 copy number variations also and make another simulated sequence. 
```bash
$ perl simuG.pl -refseq Simulation_SNP_4000_INDEL_4000.simseq.genome.fa -cnv_count 4 -cnv_min_size 20000 -prefix Simulation_SNP_4000_INDEL_4000_CNV_4 -seed 512345678900

$ ls -1tr Simulation_SNP_4000_INDEL_4000_CNV_4*
Simulation_SNP_4000_INDEL_4000_CNV_4.refseq2simseq.CNV.vcf
Simulation_SNP_4000_INDEL_4000_CNV_4.refseq2simseq.map.txt
Simulation_SNP_4000_INDEL_4000_CNV_4.simseq.genome.fa

$ sed -i '/>/ s/\(.*\)/\1_CNV_4/' Simulation_SNP_4000_INDEL_4000_CNV_4.simseq.genome.fa

$ less Simulation_SNP_4000_INDEL_4000_CNV_4.simseq.genome.fa
>NC_017518.1_SNP_4000_INDEL_4000_CNV_4
TTCGGCTTAAACCTTATCCATATCCAAACGCATAACCGTAACCCATTCACCGTTATGGAAATGTCGCCCGACAACCGCCCAGCCGAATGATTCATAAAATATTTGCACATCAGGCGTATAAAGATACAAGAACTTTATCCCCAGCGAACGCGCTGCGCCTATGCAGTGGGCGACCAGCCTCCTGCCAATGCCTTGTCCGCGATATTCAGGTAAAACAAAGACATCAC

```

### 4. Concatenate all files and make an index
We'll concatenate all these files into one and make an index too
```bash
$ cat GCF_000191525.1_ASM19152v1_genomic.fna Simulation_SNP_5000.simseq.genome.fa Simulation_INDEL_5000.simseq.genome.fa Simulation_SNP_4000_INDEL_4000.simseq.genome.fa Simulation_SNP_4000_INDEL_4000_INV_4.simseq.genome.fa Simulation_SNP_4000_INDEL_4000_CNV_4.simseq.genome.fa > ASM19152v1_pgsim.fa
$ module load SAMtools
$ samtools faidx ASM19152v1_pgsim.fa 
$ cat ASM19152v1_pgsim.fa.fai 
NC_017518.1     2248966 64      80      81
NC_017518.1_SNP_5000    2248966 2277165 2248966 2248967
NC_017518.1_INDEL_5000  2249048 4526156 2249048 2249049
NC_017518.1_SNP_4000_INDEL_4000 2242147 6775238 2242147 2242148
NC_017518.1_SNP_4000_INDEL_4000_INV_4   2242147 9017425 2242147 2242148
NC_017518.1_SNP_4000_INDEL_4000_CNV_4   2415498 11259612        2415498 2415499
```

