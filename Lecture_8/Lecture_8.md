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

*Remember to connect to an open node (i.e. pinnacle-l6)

```
#It is not available as module
#Let's activate our conda environment

module load python/anaconda-3.8
source /share/apps/bin/conda-3.8.sh

conda activate fastqc

conda install -c bioconda seqtk

```

Now we are going to use Seqtk to manipulate our fastq file:

```
#Subsample 10 reads
seqtk sample -s100 Sacc.R1_paired.trimmed.fastq.gz 10 > Sacc_sub1.fq

#Convert reads to fasta
 seqtk seq -aQ64 -q20 Sacc_sub1.fq > Sacc_sub1.fa
```

With this ready now we can use blast, but before there are some key arguments to use for analysis.
You can see a full list typing `blastn -help` to print a detailed list.

### BLAST: some key arguments
* -query - query file name (required)
* -db - database file name (require)
* -evalue - set the evalue cutoff
* -max_target_seqs - max number of hit seqs to show
* -num_alignments - max number of alignments to show
* -num_threads - number of threads (parallel processing to run, 8 will be
faster than 2)
* -outfmt - specify a simpler format than the text format, try ‘-outfmt 6’ for
tabular format
* -subject - instead of doing a DB search, search for alignments between
query sequence and 1 to many subject sequences. Useful when want to just
see the alignment of 2 sequences already picked out from other analyses

```
blastn -query Sacc_sub1.fa -db Sb

```

We can modify the types of outputs of the blast, one of the most useful output is a tabulated format, 
by using a `-outfmt 6` or `-outfmt 7`, the standard output will include:

```
 
   1.  qseqid      query or source (e.g., gene) sequence id
   2.  sseqid      subject  or target (e.g., reference genome) sequence id
   3.  pident      percentage of identical matches
   4.  length      alignment length (sequence overlap)
   5.  mismatch    number of mismatches
   6.  gapopen     number of gap openings
   7.  qstart      start of alignment in query
   8.  qend        end of alignment in query
   9.  sstart      start of alignment in subject
  10.  send        end of alignment in subject
  11.  evalue      expect value
  12.  bitscore    bit score

```

But it can be redefined to `-outfmt "6 qseqid sseqid evalue"`, limiting the number of columns.

How do submit a script `run_blast.sh`, let's put it together:

```
#!/bin/bash
#SBATCH --job-name=blast_jarojas
#SBATCH --output=blas_job1.slurm
#SBATCH --nodes=1
#SBATCH --tasks-per-node=4
$SBATCH --time=00:05:00
#SBATCH --partition comp06

module blast/2.12.0+bin

blastn -query Sacc_sub1.fa -db Sb -outfmt 6 -out reads-vs-Sbayanus.blastn.tab -num_threads 4

```

Now you can submit the job and check on it:

* `$ sbatch run_blast.sh` Actual submission
* `$ squeue -u $USER`  Check the submitted job
