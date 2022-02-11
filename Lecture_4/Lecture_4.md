---
title: "Lecture 4"
author: "Alejandro Rojas"
date: "2/11/2022"
output: 
  html_document: 
    keep_md: yes
---



## Installing Git

__MacOs__
A different option to git on mac apart from Xcode - command line, you can also install it via a binary installer. A macOS Git installer is maintained and available for download at the Git website, at https://git-scm.com/download/mac.

__Windows__

There are also a few ways to install Git on Windows. The most official build is available for download on the Git website. Just go to https://git-scm.com/download/win and the download will start automatically. Note that this is a project called Git for Windows, which is separate from Git itself; for more information on it, go to https://gitforwindows.org.

After git is installed, please check:
```
git --version
```

## Git - How it works?

We’ll start by exploring how version control can be used to keep track of what one person did and when. Even if you aren’t collaborating with other people, automated version control is much better than this situation:

"Piled Higher and Deeper" by Jorge Cham, http://www.phdcomics.com

![][id1]

Version control systems start with a base version of the document and then record changes you make each step of the way. You can think of it as a recording of your progress: you can rewind to start at the base document and play back each change you made, eventually arriving at your more recent version.

Once you think of changes as separate from the document itself, you can then think about “playing back” different sets of changes on the base document, ultimately resulting in different versions of that document. For example, two users can make independent sets of changes on the same document. Unless multiple users make changes to the same section of the document - a conflict - you can incorporate two sets of changes into the same base document.

![][id2]



[id1]: Images/phd101212s.png
[id2]: Images/Version_control.png
