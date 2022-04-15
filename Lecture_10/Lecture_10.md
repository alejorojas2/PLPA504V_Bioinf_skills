---
title: "Lecture_10"
author: "Juanita Gil"
date: "4/15/2022"
output: 
  html_document: 
    keep_md: yes
---


# Viewing Alignments

Now that we have aligned reads to a reference genome after quality review and trimming, let’s look at the alignment files. We are going to use the Interactive Genome Viewer (IGV), a tool that allows us to view alignments in an interactive way.

First, go to [pinnacle portal](https://pinnacle-portal.uark.edu/pun/sys/) dashboard and click on Interactive Apps. Select IGV from the pull-down menu. You should see the following:


![][id1]


In the different options, let's change the parameters:
- Hours = 6
- Cores to 4
- Queue to comp06 

Then click Launch!

Your window should have changed to this once the session starts:


![][id2]



Now you click in the blue button that  says __"Launch IGV"__
Then you will see the screen with IGV running, and it is GUI (Graphical User Interface), where we can start evaluating the files derived from the mapping.


![][id3]


## Uploading the files to IGV

### Upload reference genome:

Go to Genomes -> Load Genome from File.. and select your reference genome (.fna file) from the data folder. Click Open.


### Upload alignment file:

Go to Files -> Load from File.. and select your sorted bam file from the results/bam folder. Click Open.

### Upload annotation file:
Go to Files -> Load from File.. and select the gff file from your data folder. Click Open.

Once you have all the files loaded, to view the alignments, select one of the chromosomes from the pull-down menu. Use the `+` button to zoom in and view the reads.


![][id4]


There are multiple panels and each of them represent different tracks:

* Track 1: The upper panel with vertical gray bars shows depth of coverage per variant site. Pointing the mouse arrow on each bar will open a pop-up box that shows the actual depth of coverage (number of reads) at the specific site, and the number and percentage of reads with a specific allele.

* Track 2: This track shows the aligned reads. The direction at which the arrows point represents Forward (=>) and Reverse (<=) reads.

* Track 3: This track shows the genes annotated in the reference genome.


### Interpreting the colors in IGV

* __SNPs__ are shown in different colors according to the base.

* __Small insertions__ are indicated by a purple ‘I’. When the mouse arrow is pointed on the insertion, a pop-up box shows the inserted bases. 

* When the arrow is pointed at a specific read, a pop-up box appears with details on the specific read, e.g., read name, alignment start, mapping quality, base phred quality, insert size, etc. 

* At high resolution, the __specific base or variant site__ that falls at the center of the IGV window is marked by two vertical dashed lines. When a variant site has an allele that is different from reference in more than 25% of reads, the alternate allele is indicated in a different color. 

* A __deletion__ is indicated by a black horizontal bar within gaps in a read. 

* An __insert size__ that is smaller than expected is represented as a blue bar.

* Reads that are colored red have larger than expected inferred sizes.


[id1]: Images/Step1.png
[id2]: Images/Step2.png
[id3]: Images/Step3.png
[id4]: Images/Step4.png
