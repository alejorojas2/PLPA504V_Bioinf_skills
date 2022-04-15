---
title: "Lecture_10"
author: "Alejandro Rojas"
date: "4/15/2022"
output: 
  html_document: 
    keep_md: yes
---

Let's start an interactive session:

```
srun --nodes=1 --ntasks-per-node=1  --cpus-per-task=8 --partition cloud72 --time=6:00:00 --pty /bin/bash
```

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



[id1]: Images/Step1.png
[id2]: Images/Step2.png
[id3]: Images/Step3.png
[id4]: Images/Step4.png
