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
mkdir -p data results/sam results/bam results/bcf results/vcf
```

Also, download our reference data:
```
git clone git@github.com:alejorojas2/PLPA504V_Bioinf_example_data.git
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
bwa mem ref_genome.fasta input_file_R1.fastq input_file_R2.fastq > results/sam/Sacc_bayanu.S1.aligned.sam

```

The SAM file, is a tab-delimited text file that contains information for each individual read and its alignment to the genome. While we do not have time to go into detail about the features of the SAM format, the paper by [Heng Li et al.](http://bioinformatics.oxfordjournals.org/content/25/16/2078.full) provides a lot more detail on the specification.

The compressed binary version of SAM is called a BAM file. We use this version to reduce size and to allow for indexing, which enables efficient random access of the data contained within the file.

The file begins with a header, which is optional. The header is used to describe the source of data, reference sequence, method of alignment, etc., this will change depending on the aligner being used. Following the header is the alignment section. Each line that follows corresponds to alignment information for a single read. Each alignment line has 11 mandatory fields for essential mapping information and a variable number of other fields for aligner specific information. An example entry from a SAM file is displayed below with the different fields highlighted.

![][id1]

We will convert the SAM file to BAM format using the samtools program with the view command and tell this command that the input is in SAM format (-S) and to output BAM format (-b):

```

module load samtools

# Sbayanu
samtools view -S -b sam/Sacc_bayanu.S1.aligned.sam > bam/Sacc_bayanu.S1.aligned.bam

# Scerevisae
samtools view -S -b sam/Sacc_cerevisae.S1.aligned.sam > bam/Sacc_cerevisae.S1.aligned.bam

```

Next we sort the BAM file using the sort command from samtools. `-o` tells the command where to write the output.

```
# Sbayanu
samtools sort -o bam/Sacc_bayanu.S1.aligned.sorted.bam bam/Sacc_bayanu.S1.aligned.bam

# Scerevisae
samtools sort -o bam/Sacc_cerevisae.S1.aligned.sorted.bam bam/Sacc_cerevisae.S1.aligned.bam
```

You can use samtools to learn more about this bam file as well.

```
# Sbayanu
samtools flagstat bam/Sacc_bayanu.S1.aligned.sorted.bam

```
This is the output for _S. bayanu_:

```
132126 + 0 in total (QC-passed reads + QC-failed reads)
0 + 0 secondary
6 + 0 supplementary
0 + 0 duplicates
34396 + 0 mapped (26.03% : N/A)
132120 + 0 paired in sequencing
66060 + 0 read1
66060 + 0 read2
30700 + 0 properly paired (23.24% : N/A)
31144 + 0 with itself and mate mapped
3246 + 0 singletons (2.46% : N/A)
100 + 0 with mate mapped to a different chr
34 + 0 with mate mapped to a different chr (mapQ>=5)
```

Now, let's repeat it for _S. cerevisae_:

```
# Scerevisae
samtools flagstat bam/Sacc_cerevisae.S1.aligned.sorted.bam

#Output
132198 + 0 in total (QC-passed reads + QC-failed reads)
0 + 0 secondary
78 + 0 supplementary
0 + 0 duplicates
126264 + 0 mapped (95.51% : N/A)
132120 + 0 paired in sequencing
66060 + 0 read1
66060 + 0 read2
124916 + 0 properly paired (94.55% : N/A)
126042 + 0 with itself and mate mapped
144 + 0 singletons (0.11% : N/A)
122 + 0 with mate mapped to a different chr
60 + 0 with mate mapped to a different chr (mapQ>=5)

```

This is our pipeline for the analysis so far:

![][id2]



# Variant calling

First we need to get bcftools in the server:

```
#It is not available as module
#Let's activate our conda environment

module load python/anaconda-3.8
source /share/apps/bin/conda-3.8.sh

conda activate fastqc

conda install -c bioconda bcftools
```

Do the first pass on variant calling by counting read coverage with bcftools. We will use the command mpileup. The flag -O b tells bcftools to generate a bcf format output file, -o specifies where to write the output file, and -f flags the path to the reference genome:


```
bcftools mpileup -O b -o bcf/Sacc_bayanu.bcf \
          -f ~/PLPA504V_Bioinf_example_data/Saccharomyces_bayanus//Sbayanus_ASM1943126v1_genomic.fna \
          bam/Sacc_bayanu.S1.aligned.sorted.bam

```

Identify SNVs using bcftools call. We have to specify ploidy with the flag `--ploidy`, which is one for the haploid _S. cerevisae_. 

`-m` allows for multiallelic and rare-variant calling, `-v` tells the program to output variant sites only (not every site in the genome), and `-o` specifies where to write the output file:

```
bcftools call --ploidy 1 -m -v -o results/vcf/Sacc_bayanu_variants.vcf results/bcf/Sacc_bayanu.bcf 

```
[id1]: Images/sam_bam.png
[id2]: Images/pipeline.png


