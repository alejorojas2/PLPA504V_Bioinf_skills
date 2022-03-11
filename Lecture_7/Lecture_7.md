---
title: "Lecture_7"
author: "Alejandro Rojas"
date: "3/11/2022"
output: 
  html_document: 
    keep_md: yes
---


# QC for reads

Let's start an interactive session:

```
srun --nodes=1 --ntasks-per-node=1  --cpus-per-task=8 --partition comp06 --time=6:00:00 --pty /bin/bash
```

Last class, we were talking about [FASTQ format](https://en.wikipedia.org/wiki/FASTQ_format).

In FASTQ files, quality values are encoded by symbols and each of those symbols represent a number between 0 - 40.
Quality scores are provided per base, in a log-10 scaled format. The likelihood of a base call being erroneous is 10 to the power of (-Q/10) - so for a Q value of 3, you get an error rate of 1/2, for a Q value of 30, you get an error rate of 0.001, and for a Q value of 40, you get an error rate of .0001.

We already checked the qualities for our libraries, and we did using interactive session:

```
#First load the modules
module load java
module load fastqc

# In the folder containing the fastq files
fastqc *.fastq.gz
```

Or using a submission script, in our case we called `fastqc_job.sh` and here is what that file contains:
```
#!/bin/bash
#SBATCH --job-name=fastqc_jarojas
#SBATCH --output=fastqc_job1.slurm
#SBATCH --nodes=1
#SBATCH --tasks-per-node=8
$SBATCH --time=00:05:00
#SBATCH --partition comp06

module purge
module load java
module load fastqc

cd ~/Lecture_6
fastqc *.fastq.gz
```

what effect will losing real data have vs including spurious data, given the biological goal and the tools you have to reach it?

![][id3]

For that reason we have to do a QC check for our data.  Basically, trimming tends to decrease sensitivity, while increasing specificity.


## Trimming low quality reads

There are different options to approach this part of the analysis:
* Trimmomatic
* FASTX-Toolkit
* cutadapt

We will use a program called [Trimmomatic](http://www.usadellab.org/cms/?page=trimmomatic) to filter poor quality reads and trim poor quality bases from our samples.

```
#It is not available as module
#Let's activate our conda environment

module load python/anaconda-3.8
source /share/apps/bin/conda-3.8.sh

conda activate fastqc

# Now let's intall trimmomatic

conda install -c bioconda trimmomatic
```

Let's check `trimmomatic`:

```
trimmomatic
```

This output shows us that we must first specify whether we have paired end (PE) or single end (SE) reads. Next, we specify what flag we would like to run. For example, you can specify threads to indicate the number of processors on your computer that you want Trimmomatic to use. In most cases using multiple threads (processors) can help to run the trimming faster. These flags are not necessary, but they can give you more control over the command. The flags are followed by positional arguments, meaning the order in which you specify them is important. In paired end mode, Trimmomatic expects the two input files, and then the names of the output files. These files are described below. While, in single end mode, Trimmomatic will expect 1 file as input, after which you can enter the optional settings and lastly the name of the output file.

![][id4]

The last thing trimmomatic expects to see is the trimming parameters:

![][id5]

Let's try this in our files:

```
trimmomatic PE -threads 4 Sacc.R1.fastq Sacc.R2.fastq  \
              Sacc.R1.trimmed.fastq Sacc.R1un.trimmed.fastq \
              Sacc.R2.trimmed.fastq Sacc.R2un.trimmed.fastq \
              ILLUMINACLIP:SRR_adapters.fa SLIDINGWINDOW:4:20
```

[id3]: Images/Illumina.png
[id4]: Images/Parameters_1.png
[id5]: Images/Parameters_2.png
