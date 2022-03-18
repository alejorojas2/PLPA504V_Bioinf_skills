---
title: "Lecture_8"
author: "Alejandro Rojas"
date: "3/18/2022"
output: 
  html_document: 
    keep_md: yes
---

Let's start an interactive session:

```
srun --nodes=1 --ntasks-per-node=1  --cpus-per-task=8 --partition cloud72 --time=6:00:00 --pty /bin/bash
```

# Sequence search tools - BLAST
* BLAST is by far the most taught tool in Bioinformatics. I am not going to rehash this completely in this class.
* See NCBI’s [Introduction to BLAST](https://www.ncbi.nlm.nih.gov/books/NBK1762/)
* One of 7 Million pages by Googling “blast introduction tutorial”

# BLAST on HPCC
There are multiple flavors of BLAST (implementations). Focus on the latest version from NCBI (2.12.0+). Default on the cluster is blast/2.12.0+bin

```
# Checking blast modules available
module spider blast

# Loading the most recent blast
module blast/2.12.0+bin
```

First let's clone the data into your home folder from a github: https://github.com/alejorojas2/PLPA504V_Bioinf_example_data

```
git clone git@github.com:alejorojas2/PLPA504V_Bioinf_example_data.git
```

We will make links to a file which is ORFs from one yeast species _Saccharomyces bayanu_

```
#this is similar to an alias or quick access
ln -s ~/PLPA504V_Bioinf_example_data/Saccharomyces_bayanus/
```

Now we have some files, set them up for running BLAST. Our question is, what
ORFs are similar at the DNA level between these two species.

```
#Let's extract the fasta file first
gunzip Sbayanus_ASM1943126v1_genomic.fna.gz

#Let's make a database
makeblastdb -dbtype nucl -in Sbayanus_ASM1943126v1_genomic.fna
```

Now, we need to create a reduced file to check our sequences.  Let's install [seqtk](https://docs.csc.fi/apps/seqtk/) for manipulating fasta and fastq files:

```
#It is not available as module
#Let's activate our conda environment

module load python/anaconda-3.8
source /share/apps/bin/conda-3.8.sh

conda install -c bioconda seqtk

```

Now we are going to use Seqtk to manipulate our fastq file:

```
#Subsample 10 reads
seqtk sample -s100 Sacc.R1_paired.trimmed.fastq.gz 10 > Sacc_sub1.fq

#Convert reads to fasta
 seqtk seq -aQ64 -q20 Sacc_sub1.fq > Sacc_sub1.fa
```

With this ready now we can use blast:

```
blastn -query Sacc_sub1.fa -db Sbayanus_ASM1943126v1_genomic.fna

```