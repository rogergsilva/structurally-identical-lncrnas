# Exploring the Hidden Hot World of Long Non-coding RNAs in Thermophilic Fungus Using a Robust Computational Pipeline.

<p align="center">
  <img src="https://github.com/rogergsilva/structurally-identical-lncrnas/blob/main/images/image16.jpg" width="250" alt="Structurally identical tanscripts">
</p>

- [Abstract](#Abstract)
- [Transcriptome pipeline](#Transcriptome-pipeline)
  - [Requirements](#Requirements)
    - [Step 1: download and install miniconda](#step-1-download-and-install-miniconda)
    - [Step 2: setting up bioconda channels](#step-2-setup-bioconda-channels)
    - [Step 3: download and build fastp](#step-3-download-and-build-fastp)
- [Bugs](#Bugs)
- [Citation](#Citation)

# Abstract

Long noncoding RNAs (lncRNAs) are versatile RNA molecules recently identified as key regulators of gene expression in response to environmental stress. Our primary focus in this study was to develop a robust computational pipeline for identifying structurally identical lncRNAs across replicates from publicly available bulk RNA-seq datasets. To demonstrate the effectiveness of the pipeline, we utilized the transcriptome of the thermophilic fungus Thermothelomyces thermophilus and assessed the expression pattern of lncRNA in conjunction with Heat Shock Proteins (HSP), a well-known family critical for the cell's response to elevated temperatures. Our results demonstrate that the identification of structurally identical transcripts among replicates in this thermophilic fungus ensures the reliability and accuracy of RNA studies, contributing to the validity of biological interpretations. Furthermore, the majority of lncRNAs exhibited a distinct expression pattern compared to HSPs. We hope this investigation contributes to advancing our comprehension of the biological mechanisms involving lncRNAs in thermophilic fungi.

# Transcriptome pipeline

## Requirements
* Internet access;
* Python installed (> 3.10) - https://www.python.org/

### Step 1: Download and install miniconda
Please, see https://docs.anaconda.com/free/miniconda/index.html

### Step 2: Setup bioconda channels
After having installed miniconda, you will need to set up the bioconda channels.
Follow the Usage steps on https://bioconda.github.io/

### Step 3: Create a folder for your project in your computer
Create a folder in your computer where you will carry out your in-silico experiments.
For example, if your project will name PRJ123, you will execute the following command:
Linux operation system
```shell
mkdir PRJ123
```
### Step 4: Download the pipeline script file
Download the file pipeline.py from the Github and save it in the folder you have created. 
For example, save it in the folder PRJ123.

### Step 5: Testing the pipeline script file
Open a command line terminal, change the path to the folder PRJ123 and execute the following command:
```shell
cd PRJ123
python pipeline.py --about
```
If you see something like: Transcriptome assembly pipeline - version X

### Step 6: Creating required folders
On the command line terminal, change the path to the folder PRJ123 and execute the following command:
```shell
cd PRJ123
python pipeline.py --about
```

# Bugs
Do you have an issue? Please file it on :https://github.com/rogergsilva/structurally-identical-lncrnas/issues

# Citation
