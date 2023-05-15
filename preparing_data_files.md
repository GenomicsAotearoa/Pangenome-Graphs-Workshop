# Preparing required data files with simuG tool
[simuG](https://github.com/yjx1217/simuG) is a light-weighted tool which can be used to simulate sequences with pre-defined number of random genomic variants against a reference sequence (Please refer the [Git Page](https://github.com/yjx1217/simuG) for more details). 

### 1. Setting up the tool 
Please refer the installation instruction in the [Git Page](https://github.com/yjx1217/simuG).
```bash
git clone https://github.com/yjx1217/simuG.git
cd simuG
perl simuG.pl -h
perl vcf2model.pl -h
```

### 2. Data source
You can download the reference sequence (Neisseria gonorrhoeae FA 1090 genome assembly) from the site [Genome assembly ASM684v1](https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_000006845.1/) and get the file _\ncbi_dataset\data\GCA_000006845.1\GCA_000006845.1_ASM684v1_genomic.fna_ from he zip file ncbi_dataset.zip. Rename it as "ref.fa" to refere easily in next steps. 

In Unix environment you can use `curl`. 

```bash
$ curl -OJX GET "https://api.ncbi.nlm.nih.gov/datasets/v2alpha/genome/accession/GCF_000006845.1/download?include_annotation_type=GENOME_FASTA,GENOME_GFF,RNA_FASTA,CDS_FASTA,PROT_FASTA,SEQUENCE_REPORT&filename=GCF_000006845.1.zip" -H "Accept: application/zip"
  % Total    % Received % Xferd  Average Speed   Time    Time     Time  Current
                                 Dload  Upload   Total   Spent    Left  Speed
100 1807k    0 1807k    0     0   294k      0 --:--:--  0:00:06 --:--:--  421k

$ ls
GCF_000006845.1.zip

$ unzip GCF_000006845.1.zip 
Archive:  GCF_000006845.1.zip
  inflating: README.md               
  inflating: ncbi_dataset/data/assembly_data_report.jsonl  
  inflating: ncbi_dataset/data/GCF_000006845.1/GCF_000006845.1_ASM684v1_genomic.fna  
  inflating: ncbi_dataset/data/GCF_000006845.1/genomic.gff  
  inflating: ncbi_dataset/data/GCF_000006845.1/cds_from_genomic.fna  
  inflating: ncbi_dataset/data/GCF_000006845.1/protein.faa  
  inflating: ncbi_dataset/data/GCF_000006845.1/sequence_report.jsonl  
  inflating: ncbi_dataset/data/dataset_catalog.json  

$ cp ncbi_dataset/data/GCF_000006845.1/GCF_000006845.1_ASM684v1_genomic.fna ref.fa

$ head ref.fa 
>NC_002946.2 Neisseria gonorrhoeae FA 1090, complete sequence
ATAAATTTTTGCACGGGTTGTGGATAAAATATCGGCGAGTCGGTATAATCGGTTCGCTGCGTTTTGAACCGACGCGTATT
CAACAGATTTGTTTTCTTTTTGAAAATATTATATTTTCTTTGTTTTCGATTTCATTTTTACCGATTCGAGCCTATCGCAT
GACATTAGCAGAGTTTTGGCCGCTGTGCCTCCGCCGTCTTCACGATATGTTGCCTCACGGGCAGTTTGCGCAATGGATTG
CGCCCCTTACGGTTGGTGAGGAGGGTGGCGTATGGGTGGTGTACGGCAAGAACCAGTTTGCCTGCAATATGCTCAAGAGC
CAGTTTGCCGGAAAAATAGAAGCGGTAAGGGAAGAGTTGGCTGCCGGCCGTCCCGCCTTTGTATTCAAACCGGGAGAAGG
CGTGCGTTATGAGATGGCGGCGGTTGAAGGTGCTGTCGAACCTGCCGAGCCGTCCTTGCACGCGGGGTCGGAGGAGATGC
CCGTGCAGGAGGTTCTGTTGGACGAGCTGCCGTCTGAAAAGCCTGTCAAACCCGCTGCGTCGAAAACGGCGGCGGATATT
TTGGCGGAACGTATGAAAAACCTGCCGCACGAGCCGCGTCAGGCTGCCGGGCCTGCTTCCCGGCCGGAATCGGCGGCAGT
TGCCAAAGCGCGGACGGATGCGCAGCGTGATGCGGAAGAAGCGCGTTACGAACAAACCAACCTGTCTCCGGATTACACGT
```

### 3. Simulating sequences with Known SNPs and INDELs counts using simuG
#### i. SNP=3000 and INDEL=300
```
$ perl simuG.pl -refseq ref.fa -snp_count 3000 -indel_count 300 -prefix sim3k -seed 123456789000

$ ls -1sh sim3k*
256K sim3k.refseq2simseq.INDEL.vcf
512K sim3k.refseq2simseq.map.txt
512K sim3k.refseq2simseq.SNP.vcf
2.3M sim3k.simseq.genome.fa
```

#### ii. SNP=4000 and INDEL=400
```bash
$ perl simuG.pl -refseq ref.fa -snp_count 4000 -indel_count 400 -prefix sim4k -seed 123456789000

$ ls -1sh sim4k*
256K sim4k.refseq2simseq.INDEL.vcf
512K sim4k.refseq2simseq.map.txt
768K sim4k.refseq2simseq.SNP.vcf
2.3M sim4k.simseq.genome.fa
```

#### iii. SNP=5000 and INDEL=500
```bash
$ perl simuG.pl -refseq ref.fa -snp_count 5000 -indel_count 500 -prefix sim5k -seed 123456789000

$ ls -1sh sim5k*
256K sim5k.refseq2simseq.INDEL.vcf
512K sim5k.refseq2simseq.map.txt
768K sim5k.refseq2simseq.SNP.vcf
2.3M sim5k.simseq.genome.fa
```

### 4. Renaming file headers for clarity
Reference file has long header name and simulated files have same header name. We have to rename them for clarity. `sed` can be used for that. 
```
$ head ref.fa 
>NC_002946.2 Neisseria gonorrhoeae FA 1090, complete sequence
ATAAATTTTTGCACGGGTTGTGGATAAAATATCGGCGAGTCGGTATAATCGGTTCGCTGCGTTTTGAACCGACGCGTATT
CAACAGATTTGTTTTCTTTTTGAAAATATTATATTTTCTTTGTTTTCGATTTCATTTTTACCGATTCGAGCCTATCGCAT
GACATTAGCAGAGTTTTGGCCGCTGTGCCTCCGCCGTCTTCACGATATGTTGCCTCACGGGCAGTTTGCGCAATGGATTG
CGCCCCTTACGGTTGGTGAGGAGGGTGGCGTATGGGTGGTGTACGGCAAGAACCAGTTTGCCTGCAATATGCTCAAGAGC
CAGTTTGCCGGAAAAATAGAAGCGGTAAGGGAAGAGTTGGCTGCCGGCCGTCCCGCCTTTGTATTCAAACCGGGAGAAGG
CGTGCGTTATGAGATGGCGGCGGTTGAAGGTGCTGTCGAACCTGCCGAGCCGTCCTTGCACGCGGGGTCGGAGGAGATGC
CCGTGCAGGAGGTTCTGTTGGACGAGCTGCCGTCTGAAAAGCCTGTCAAACCCGCTGCGTCGAAAACGGCGGCGGATATT
TTGGCGGAACGTATGAAAAACCTGCCGCACGAGCCGCGTCAGGCTGCCGGGCCTGCTTCCCGGCCGGAATCGGCGGCAGT
TGCCAAAGCGCGGACGGATGCGCAGCGTGATGCGGAAGAAGCGCGTTACGAACAAACCAACCTGTCTCCGGATTACACGT

$ sed -i '/>/ s/.*/>ref/g' ref.fa

$ >ref
ATAAATTTTTGCACGGGTTGTGGATAAAATATCGGCGAGTCGGTATAATCGGTTCGCTGCGTTTTGAACCGACGCGTATT
CAACAGATTTGTTTTCTTTTTGAAAATATTATATTTTCTTTGTTTTCGATTTCATTTTTACCGATTCGAGCCTATCGCAT
GACATTAGCAGAGTTTTGGCCGCTGTGCCTCCGCCGTCTTCACGATATGTTGCCTCACGGGCAGTTTGCGCAATGGATTG
CGCCCCTTACGGTTGGTGAGGAGGGTGGCGTATGGGTGGTGTACGGCAAGAACCAGTTTGCCTGCAATATGCTCAAGAGC
CAGTTTGCCGGAAAAATAGAAGCGGTAAGGGAAGAGTTGGCTGCCGGCCGTCCCGCCTTTGTATTCAAACCGGGAGAAGG
CGTGCGTTATGAGATGGCGGCGGTTGAAGGTGCTGTCGAACCTGCCGAGCCGTCCTTGCACGCGGGGTCGGAGGAGATGC
CCGTGCAGGAGGTTCTGTTGGACGAGCTGCCGTCTGAAAAGCCTGTCAAACCCGCTGCGTCGAAAACGGCGGCGGATATT
TTGGCGGAACGTATGAAAAACCTGCCGCACGAGCCGCGTCAGGCTGCCGGGCCTGCTTCCCGGCCGGAATCGGCGGCAGT
TGCCAAAGCGCGGACGGATGCGCAGCGTGATGCGGAAGAAGCGCGTTACGAACAAACCAACCTGTCTCCGGATTACACGT
```

Renaming simulated files
```bash
$ sed -i '/>/ s/.*/>sim3k/g' sim3k.simseq.genome.fa
$ sed -i '/>/ s/.*/>sim4k/g' sim4k.simseq.genome.fa
$ sed -i '/>/ s/.*/>sim5k/g' sim5k.simseq.genome.fa
```

### 5. Concatenate all file and make an index
We'll concatenate these files into one and make an index too
```bash
$ cat ref.fa sim3k.simseq.genome.fa sim4k.simseq.genome.fa sim5k.simseq.genome.fa > 4Sim.fa
$ module load SAMtools
$ samtools faidx 4Sim.fa 
$ cat 4Sim.fa.fai 
ref     2153922 5       80      81
sim3k   2153907 2180859 2153907 2153908
sim4k   2153912 4334774 2153912 2153913
sim5k   2153883 6488694 2153883 2153884
```

