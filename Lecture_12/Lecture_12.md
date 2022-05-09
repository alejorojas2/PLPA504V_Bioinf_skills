---
title: "Lecture_12"
author: "Alejandro Rojas"
date: "5/9/2022"
output: 
  html_document: 
    keep_md: yes
    theme: readable
---

# Downloading data from SRA

Most scientific journals require scientists to make their sequencing data publicly available. This way, other researchers in the world can download the raw data and re-analyze it for their own purposes. 

How do we access the data? Raw sequencing data comes in huge files that are often multiple gigabytes in size per sample. If you are a researcher with little bioinformatics experience, the finding and downloading the data can be somewhat complicated. This guide explains how to:

1. Navigate through SRA to find raw sequencing data.
2. ownload and convert SRA files to FASTQ files using the NCBI's SRA toolkit.

## Finding a dataset in SRA

Let's say you are reading a paper in a journal and see an interesting RNA-seq experiment. You decide that you want to sift through the data for your own genes of interest. The first step is finding the Bioproject or Sample accession number corresponding to the dataset.

![][id1]

## Downloading a SRA dataset

In order to download the SRA files onto your machine, we use the NCBI's SRA toolkit, which lets us use the command line to download a specified SRA run. If you are using the UARK HPCC, you need to load the modules.

```
# Loading modules

module load java
module load sra-toolkit

# Module check (Print version of Prefetch, parto of sra-toolkit)
prefetch -V
```

You can read more about SRA toolkit here: https://www.ncbi.nlm.nih.gov/books/NBK242621/ and at their github repo: https://github.com/ncbi/sra-tools.

The toolkit works by first using the prefetch command to download the SRA file associated with the specified SRA run ID.

For example, to download the SRA file for HET_CD4_1 (SRA Run identifier: SRR2121685), the command would be: 

```
prefetch SRR2121685
```

__However, right now it doesn't work in server__, but the expected outcome is a single file `SRR2121685.sra`, if it is working you should see:

```
2020-02-06T21:54:29 prefetch.2.8.2: 1) Downloading 'SRR2121685'...
2020-02-06T21:54:29 prefetch.2.8.2:  Downloading via https...
2020-02-06T21:57:32 prefetch.2.8.2: 1) 'SRR2121685' was downloaded successfully

```

The file SRR2121685.sra should be downloaded into your home directory at `~/ncbi/public/sra/`. You can double check by listing the files:

```
ls ~/ncbi/public/sra
SRR2121685.sra

```

After you have downloaded the SRA file, you can use the command fastq-dump to extract the contents of it into a .fastq.gz file.

```
fastq-dump --outdir fastq --gzip --split-3  ~/ncbi/public/sra/SRR2121685.sra
```

If successful, you should see the following output show up in your terminal:

```
Read 27928438 spots for /home/ericklu/ncbi/public/sra/SRR2121685.sra
Written 27928438 spots for /home/ericklu/ncbi/public/sra/SRR2121685.sra
```

We can check the folder fastq/ to make sure our files were downloaded correctly by listing files using `ls fastq`


This shoould be the output:

```
SRR2121685_pass_1.fastq.gz  SRR2121685_pass_2.fastq.gz
```

We observe that two fastq files have been extracted from SRR2121685.sra. This is because the original data was produced from paired-end sequencing, which usually has both a Read1 file and Read2 file. fastq-dump has extracted the SRA file into two files, with suffix "_1" for paired-end read 1 and "_2" for paired-end read 2. 


[id1]: Images/Image1.jpeg
[id2]: Images/Step2.png
[id3]: Images/Step3.png
