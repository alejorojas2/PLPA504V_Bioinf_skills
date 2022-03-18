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
srun --nodes=1 --ntasks-per-node=1  --cpus-per-task=8 --partition cloud72 --time=6:00:00 --pty /bin/bash
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
trimmomatic PE -threads 4 Sacc.R1.fastq.gz Sacc.R2.fastq.gz  \
              Sacc.R1_paired.trimmed.fastq.gz Sacc.R1_unpaired.trimmed.fastq.gz \
              Sacc.R2_paired.trimmed.fastq.gz Sacc.R2_unpaired.trimmed.fastq.gz \
              SLIDINGWINDOW:4:20 \
              MINLEN:50 \
              ILLUMINACLIP:adapters.fa:2:30:10 \
              
```

The command above is a single line, but since some scripts or even the monitor could cut the command, you can use backslash `\` to break the command into multiple lines.  Below I am going break down the command in parts to understand the whole process:

* `trimmomatic PE -threads 4 Sacc.R1.fastq Sacc.R2.fastq` - using trimmomatic in pair-end (PE) mode to use 4 threads and we specify the input files in order `R1` followed by `R2`
* `Sacc.R1_paired.trimmed.fastq.gz Sacc.R1_unpaired.trimmed.fastq.gz` - this part specifies two output files for reads in processed in `R1`: paired reads and unpaired reads.
* `Sacc.R2_paired.trimmed.fastq.gz Sacc.R2_unpaired.trimmed.fastq.gz` - this part specifies two output files for reads processed in `R2`: paired reads and unpaired reads.
* `SLIDINGWINDOW:4:20` - This parameter specifies the window analysis, `4` is how many bases to include in the sliding window and the __average__ quality.
* `MINLEN:50` - This part specifies the minimum read length to keep, anything below is discarded
* `ILLUMINACLIP:adapters.fa:2:30:10` - This parameter specifies the option for trimming adapters (i.e. illumina nextera), the numbers in order are: mismatches allowed;palindrome; simple (This numbers are scores and suggested values range between 7-15 for simple and 30 for palindrome)

Ideally, this should be run as a script submitted to __SLURM__, in that case you should create a bash script using __nano__ as follows: `nano trim_job.sh` (you can use any name). the contents of that file will be:

```#!/bin/bash
#SBATCH --job-name=trim_jarojas
#SBATCH --output=trim_job1.slurm
#SBATCH --nodes=1
#SBATCH --tasks-per-node=8
$SBATCH --time=00:05:00
#SBATCH --partition comp06

module load python/anaconda-3.8
source /share/apps/bin/conda-3.8.sh

conda activate fastqc

cd ~/Lecture_6

trimmomatic -threads 4 Sacc.R1.fastq.gz Sacc.R2.fastq.gz \
                       Sacc.R1_illumina.trimmed.fastq.gz Sacc.R1unpaired.trimmed.fastq.gz \
                       Sacc.R2_illumina.trimmed.fastq.gz Sacc.R2unpaired.trimmed.fastq.gz \
                       SLIDINGWINDOW:4:20 \
                       MINLEN:50 \
                       ILLUMINACLIP:adapters.fa:2:30:10

```

**Important note:** For the parameter illumina clip usually trimmomatic has folder that contains different examples of adapters used for illumina libraries:

- TruSeq2 (as used in GAII machines) 
- TruSeq3 (as used by HiSeq and MiSeq machines)
- Nextera PE

Since we installed using conda, the files for those adapters are hidden in the conda environment `~/.conda/envs/fastqc/share/trimmomatic/adapters/`, however we can make a file containing all the adapters.  Create a file named `adapters.fa` containing the following sequences:

```
>TruSeq3_IndexedAdapter
AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC
>TruSeq3_UniversalAdapter
AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTA

#TruSeq3-PE-2
>PrefixPE/1
TACACTCTTTCCCTACACGACGCTCTTCCGATCT
>PrefixPE/2
GTGACTGGAGTTCAGACGTGTGCTCTTCCGATCT
>PE1
TACACTCTTTCCCTACACGACGCTCTTCCGATCT
>PE1_rc
AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTA
>PE2
GTGACTGGAGTTCAGACGTGTGCTCTTCCGATCT
>PE2_rc
AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC

#TruSeq3-PE
>PrefixPE/1
TACACTCTTTCCCTACACGACGCTCTTCCGATCT
>PrefixPE/2
GTGACTGGAGTTCAGACGTGTGCTCTTCCGATCT

TruSeq3-SE
>TruSeq3_IndexedAdapter
AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC
>TruSeq3_UniversalAdapter
AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTA

#TruSeq2-SE
>TruSeq2_SE
AGATCGGAAGAGCTCGTATGCCGTCTTCTGCTTG
>TruSeq2_PE_f
AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT
>TruSeq2_PE_r
AGATCGGAAGAGCGGTTCAGCAGGAATGCCGAG

#TruSeq2-PE
>PrefixPE/1
AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCT
>PrefixPE/2
CAAGCAGAAGACGGCATACGAGATCGGTCTCGGCATTCCTGCTGAACCGCTCTTCCGATCT
>PCR_Primer1
AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCT
>PCR_Primer1_rc
AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATT
>PCR_Primer2
CAAGCAGAAGACGGCATACGAGATCGGTCTCGGCATTCCTGCTGAACCGCTCTTCCGATCT
>PCR_Primer2_rc
AGATCGGAAGAGCGGTTCAGCAGGAATGCCGAGACCGATCTCGTATGCCGTCTTCTGCTTG
>FlowCell1
TTTTTTTTTTAATGATACGGCGACCACCGAGATCTACAC
>FlowCell2
TTTTTTTTTTCAAGCAGAAGACGGCATACGA

#Nextera
>PrefixNX/1
AGATGTGTATAAGAGACAG
>PrefixNX/2
AGATGTGTATAAGAGACAG
>Trans1
TCGTCGGCAGCGTCAGATGTGTATAAGAGACAG
>Trans1_rc
CTGTCTCTTATACACATCTGACGCTGCCGACGA
>Trans2
GTCTCGTGGGCTCGGAGATGTGTATAAGAGACAG
>Trans2_rc
CTGTCTCTTATACACATCTCCGAGCCCACGAGAC
```


[id3]: Images/Illumina.png
[id4]: Images/Parameters_1.png
[id5]: Images/Parameters_2.png
