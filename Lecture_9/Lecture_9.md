---
title: "Lecture_9"
author: "Alejandro Rojas"
date: "4/1/2022"
output: 
  html_document: 
    keep_md: yes
---

Let's start an interactive session:

```
srun --nodes=1 --ntasks-per-node=1  --cpus-per-task=8 --partition cloud72 --time=6:00:00 --pty /bin/bash
```

Now that we have reads that are clean after the quality review and trimming we can start looking a comparing with a genome or eventually assembling a genome.

We perform read alignment or mapping to determine where in the genome our reads originated from. There are a number of tools to choose from and, while there is no gold standard, there are some tools that are better suited for particular NGS analyses. We will be using the Burrows Wheeler Aligner (BWA), which is a software package for mapping low-divergent sequences against a large reference genome.

The alignment process consists of two steps:

1. Indexing the reference genome
2. Aligning the reads to the reference genome

For this, we are going to use tow genomes to map our reads to.  These are:
- _Saccharomyces cerevisae_
- _Saccharomyces bayanu_

We need to do a bit of housekeeping before we run the actual analysis. We will also  to create directories for the results that will be generated as part of this workflow. We can do this in a single line of code, because mkdir can accept multiple new directory names as input.

```
mkdir -p results/sam results/bam results/bcf results/vcf
```

# Index the reference genomes

Our first step is to index the reference genome for use by BWA. Indexing allows the aligner to quickly find potential alignment sites for query sequences in a genome, which saves time during alignment. Indexing the reference only has to be run once. The only reason you would want to create a new index is if you are working with a different reference genome or you are using a different tool for alignment.

```
module load bwa

#Please do one at the time
bwa index ~/PLPA504V_Bioinf_example_data/Saccharomyces_bayanus/Sbayanus_ASM1943126v1_genomic.fna

#Now, Saccharomyces cerevisae
bwa index ~/PLPA504V_Bioinf_example_data/Saccharomyces_cerevisae/GCF_000146045.2_R64_genomic.fna
```

While the index is created, you will see output that looks something like this:

```
[bwa_index] Pack FASTA... 0.07 sec
[bwa_index] Construct BWT for the packed sequence...
[bwa_index] 2.90 seconds elapse.
[bwa_index] Update BWT... 0.06 sec
[bwa_index] Pack forward-only FASTA... 0.04 sec
[bwa_index] Construct SA from BWT and Occ... 1.02 sec
[main] Version: 0.7.17-r1188
[main] CMD: bwa index Sbayanus_ASM1943126v1_genomic.fna
[main] Real time: 4.161 sec; CPU: 4.100 sec
```

# Align reads to reference genome

The alignment process consists of choosing an appropriate reference genome to map our reads against and then deciding on an aligner. We will use the BWA-MEM algorithm, which is the latest and is generally recommended for high-quality queries as it is faster and more accurate.

An example of what a bwa command looks like is below. This command will not run, as we do not have the files ref_genome.fa, input_file_R1.fastq, or input_file_R2.fastq.

```
bwa mem ref_genome.fasta input_file_R1.fastq input_file_R2.fastq > output.sam

```




