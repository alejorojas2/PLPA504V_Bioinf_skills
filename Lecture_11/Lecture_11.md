---
title: "Lecture_11"
author: "Alejandro Rojas"
date: "4/29/2022"
output: 
  html_document: 
    keep_md: yes
    theme: readable
---

# Genome assembly

Here we’re going to run through some of the typical steps for taking a newly sequenced fungal isolate genome from raw fastq files through to some assemblies and discussing how we might choose which assembly to move forward with 🙂

 __disclaimer__: We are going to assemble an small and simple genome, therefore some of these methods could work for your organism of interest, but it could change.  _Please read and review the literature to decide the best path for your organism_.

_This is not an authoritative or exhaustive workflow for working with a newly sequenced genome. No such thing exists! All genomes, datasets, and goals are different, and new tools are constantly being developed. The point of this page is just to give examples of some of the things we can do. Don’t let anything here, or anywhere, constrain your science to doing only what others have done!_

## Initial preparation
Let's start an interactive session:

```
srun --nodes=1 --ntasks-per-node=1  --cpus-per-task=8 --partition cloud72 --time=6:00:00 --pty /bin/bash
```

Throughout this process we’ll be using a variety of tools. They are going to be installed using __conda__ with the specific versions used at the time this was put together.

```
# Let's create a folder for this lecture at home
mkdir Lecture_11

cd Lecture_11

```

## Data QC and Error Correction

We are going to use our reads for _S. cerevisae_ that we already processed. 

Assessing the quality of our sequence data and filtering appropriately should pretty much always be one of the first things we do with our data. A great tool we’ve already used for this is `FastQC`. We already did in a previous lecture.

Also as we’ve seen, `Trimmomatic` is a pretty flexible tool that enables us to trim up our sequences based on several quality thresholds. The summary from `FastQC` didn’t look all that terrible other than semi-low quality scores towards the end of each read, so for a first pass let’s run `Trimmomatic` to clean that up.

Now, we are going to get the data:

```
#Since we already process thre reads
ln -s ~/PLPA504V_Bioinf_example_data/Lecture_9_ready/data/reads/Sacc.R*_paired.trimmed.fastq.gz ./
```

Reads could be corrected for errors.  Even though we’ve done our best to quality trim and filter our reads, there can still be errors in there that may hurt our assembly attempts. Conceptually, read-error correction tools try to find bases that have a lot of coverage where the base is the same, but maybe in a tiny fraction of reads that base is different, and they will change that base to be what the others are.

The read-error correction we are going to use here is employed by the SPAdes assembler, and is called `BayesHammer` (initial publication [here](https://bmcgenomics.biomedcentral.com/articles/10.1186/1471-2164-14-S1-S7) if you’d like to look into more sometime).

```
# Let's activate Spades
module load spades

# Error correction
spades.py -1 Sacc.R1_paired.trimmed.fastq.gz -2 Sacc.R2_paired.trimmed.fastq.gz \
           -o spades_error_corrected_reads -t 50 -m 500 \
           --only-error-correction
```

#### Code breakdown

* `-1` and `-2` - specify the input forward and reverse reads
* `-o` - specifies the output directory where the results will be stored
* `-t` - specifies how many threads will be used (“threads” are complicated)
* `-m` - specifies the maximum amount of memory to possibly be used in gigabytes; program would quit if it needed more than that (our instances have 16 GB available to them)
* `--only-error-correction` - telling SPAdes to only run the error-correction step, and not perform an assembly

The output on the screen is:

```
System information:
  SPAdes version: 3.11.1
  Python version: 2.7.5
  OS: Linux-3.10.0-1160.2.1.el7.x86_64-x86_64-with-centos-7.9.2009-Core

Output dir: /scrfs/storage/jarojas/home/test/spades_error_corrected_reads
Mode: ONLY read error correction (without assembling)
Debug mode is turned OFF

....
....

 0:01:03.014   304M / 616M  INFO    General                 (hammer_tools.cpp          : 166)   Prepared batch 0 of 66060 reads.
  0:01:04.499   316M / 616M  INFO    General                 (hammer_tools.cpp          : 175)   Processed batch 0
  0:01:04.645   316M / 616M  INFO    General                 (hammer_tools.cpp          : 185)   Written batch 0
  0:01:04.686   260M / 616M  INFO    General                 (hammer_tools.cpp          : 270)   Correction done. Changed 2200 bases in 2075 reads.
  0:01:04.686   260M / 616M  INFO    General                 (hammer_tools.cpp          : 271)   Failed to correct 0 bases out of 9884195.
  0:01:04.687     4M / 616M  INFO    General                 (main.cpp                  : 262)   Saving corrected dataset description to /scrfs/storage/jarojas/home/test/spades_error_corrected_reads/corrected/corrected.yaml
  0:01:04.688     4M / 616M  INFO    General                 (main.cpp                  : 269)   All done. Exiting.
```

## Assembly program

Now that we have our reads quality filtered and we have run an error-correction step, we’re ready to move on to assembling them! There are [lots of assembly programs](https://en.wikipedia.org/wiki/De_novo_sequence_assemblers) out there, and once again, there is no one-size-fits-all. The reason each new assembly paper shows it doing better than all the others is not because everyone is lying about things, it’s just that the data have a lot to say about which assembler is going to work the “best” – and what’s “best” is not really a straightforward criterion to shoot for either.

### SPAdes with default settings
Let’s first run a SPAdes assembly with default settings. Note that we are providing the --only-assembler flag because we already ran the read-error correction step above, let's link those file in the main folder (Lecture_11):

```
# Let's create a symbolic link to the corrected reads
ln -s ./spades_error_corrected_reads/corrected/Sacc.R1_paired.trimmed.fastq.00.0_0.cor.fastq.gz Sacc.R1_paired.trimmed.corr.fastq.gz

ln -s ./spades_error_corrected_reads/corrected/Sacc.R2_paired.trimmed.fastq.00.0_0.cor.fastq.gz Sacc.R2_paired.trimmed.corr.fastq.gz
```

We are using the outputs from that here as the input (takes < 5 min.) to run spades:

```
spades.py -1 Sacc.R1_paired.trimmed.corr.fastq.gz -2 Sacc.R2_paired.trimmed.corr.fastq.gz -o spades-default-assembly/ -t 6 --only-assembler
```

Most of the suggestions there are defaults in the version we are using, and were run with our command above, but one (running things in --careful mode) was not. This tries to find and fix mismatches after the initial assembly is finished. Here’s how we would run one with careful mode on :

```
spades.py -1 Sacc.R1_paired.trimmed.corr.fastq.gz -2 Sacc.R2_paired.trimmed.corr.fastq.gz \
         -o spades-careful-assembly -t 6 --only-assembler \
         --careful
```


### Megahit Assembly

The [MEGAHIT documentation](https://github.com/voutcn/megahit/wiki) has an assembly-tips page that [notes here](https://github.com/voutcn/megahit/wiki/Assembly-Tips#filtering-kmin1-mer) that the default settings are largely tuned for metagenomic assembly, and that for a generic assembly (like our case here with an isolate) when there is greater than 40X coverage it is suggested to set the --min-count parameter (which deals with the frequency of kmers and filtering out reads with unique kmers) to 3.
```
#Loading megahit
module load megahit

#Simple command
megahit -1 Sacc.R1_paired.trimmed.corr.fastq.gz -2 Sacc.R2_paired.trimmed.corr.fastq.gz \
         -o megahit-min-count-3-assembly/ -t 6 --min-count 3
```

## QUAST

[QUAST](https://github.com/ablab/quast) is a really nice tool for comparing multiple assemblies, and for metagenome assemblies there is a comparable MetaQUAST. We can provide QUAST with all of our assemblies, a fasta file of our reference genome, and a .gff (general feature format) file of our reference genome that contains information about its genes. Then QUAST will compare our assemblies to the reference, and also generate several reference-independent summary statistics.

```
#Let's make the fasta file available
ln -s ~/PLPA504V_Bioinf_example_data/Lecture_9_ready/data/reference/GCF_000146045.2_R64_genomic.fna ./Scerevisae_genome.fa

# let's make out reference gff available
ln -s ~/PLPA504V_Bioinf_example_data/Lecture_9_ready/data/reference/GCF_000146045.2_R64_genomic.gff ./Scerevisae_genes.gff
```

Now we can run quast:

```
#Loading quast  
git clone git@github.com:ablab/quast.git


#Running quast
./quast/quast.py -o quast-Screvisae-out -r ./Scerevisae_genome.fa \
      -g ./Scerevisae_genes.gff -t 6 -m 1000 \
      -l "SPAdes-default, SPAdes-careful, MEGAHIT-min-count-3" \
      spades-default-assembly/contigs.fasta \
      spades-careful-assembly/contigs.fasta \
      megahit-min-count-3-assembly/final.contigs.fa

```


#### CODE BREAKDOWN

Remember the backslashes (\) don’t have anything to do with the actual command, they are just there so this is not one long line and we can copy and paste (otherwise the “enter” character that comes after them would try to execute the command on each line when pasted into our terminal.

* `-o` - setting the output directory
* `-r` reference_genome - providing a fasta file of our reference genome
* `-g` reference_genome - providing the gene coordinates of our reference (this is so QUAST can see if it finds them in our assemblies)
* `-l` "..." - within these quotes are labels for each of the input assemblies (so the output has something cleaner than file names), these need to be provided here in the same order in which their input fasta files are entered next as positional arguments
the end is 3 positional arguments for our 3 input assemblies
