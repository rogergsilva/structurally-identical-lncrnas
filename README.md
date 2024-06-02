# Exploring the Hidden Hot World of Long Non-coding RNAs in Thermophilic Fungus Using a Robust Computational Pipeline.

<p align="center">
  <img src="https://github.com/rogergsilva/structurally-identical-lncrnas/blob/main/images/image16.jpg" width="250" alt="Structurally identical tanscripts">
</p>

- [Abstract](#Abstract)
- [Transcriptome pipeline](#Transcriptome-pipeline)
  - [Requirements](#Requirements)
    - [Step 1: Download and install miniconda](#step-1-download-and-install-miniconda)
    - [Step 2: Setup bioconda channelss](#step-2-setup-bioconda-channels)
    - [Step 3: Installation and requirements](#step-3-installation-and-requeriments)
    - [Step 4: Create a folder for your project in your computer](#step-4-create-a-folder-for-your-project-in-your-computer)
    - [Step 5: Download the pipeline script file](#step-5-download-the-pipeline-script-file)
    - [Step 6: How to use it](#step-6-how-to-use-it)    
- [Input and Output](#Input-Output)
- [Example of how to use this pipeline](#Example-of-how-to-use-this-pipeline)
- [Bugs](#Bugs)
- [Citation](#Citation)

# Abstract

Long noncoding RNAs (lncRNAs) are versatile RNA molecules recently identified as key regulators of gene expression in response to environmental stress. Our primary focus in this study was to develop a robust computational pipeline for identifying structurally identical lncRNAs across replicates from publicly available bulk RNA-seq datasets. To demonstrate the effectiveness of the pipeline, we utilized the transcriptome of the thermophilic fungus Thermothelomyces thermophilus and assessed the expression pattern of lncRNA in conjunction with Heat Shock Proteins (HSP), a well-known family critical for the cell's response to elevated temperatures. Our results demonstrate that the identification of structurally identical transcripts among replicates in this thermophilic fungus ensures the reliability and accuracy of RNA studies, contributing to the validity of biological interpretations. Furthermore, the majority of lncRNAs exhibited a distinct expression pattern compared to HSPs. We hope this investigation contributes to advancing our comprehension of the biological mechanisms involving lncRNAs in thermophilic fungi.

# Identification of identical lncRNAs pipeline

## Requirements
* Internet access;
* Python installed (> 3.10) - https://www.python.org/

### Step 1: Download and install miniconda
Please, see https://docs.anaconda.com/free/miniconda/index.html

### Step 2: Setup bioconda channels
After having installed miniconda, you will need to set up the bioconda channels.
Follow the Usage steps on https://bioconda.github.io/

## Setp 3: Installation and requirements
### Requirements

The following python libraries must be installed on your machine:
```shell
  conda install statistics, biopython
```
### Step 4: Create a folder for your project in your computer
Create a folder in your computer where you will carry out your in-silico experiments.
For example, if your project will name PRJ123, you will execute the following command:
Linux operation system
```shell
mkdir PRJ123
```
### Step 5: Download the pipeline script file
Download the file identical.py from the Github and save it in the folder you have created. 
For example, save it in the folder PRJ123.

## Step 6: How to use it
Open a command line terminal, change the path to the folder PRJ123 and execute the following command:
```shell
cd PRJ123
python identical.py -h
usage: identical.py [-h] [-i IDS] -m MOLECULES -g GROUPS [-t TRACKING] [-f FASTA] [-o OUTPUT]

Identify structurally identical transcripts.

options:
  -h, --help            show this help message and exit
  -i IDS, --ids IDS     File containing only transfrag IDs that were previsouly checked - putative transcripts
  -m MOLECULES, --molecules MOLECULES
                        Which molecules will be selected to be analyzed
  -g GROUPS, --groups GROUPS
                        Group your samples in the same order they were processed by StringTie
  -t TRACKING, --tracking TRACKING
                        Tracking file from StringTie2 output to be read
  -f FASTA, --fasta FASTA
                        Fasta file from StringTie2 output to be used to save identical transcripts fasta file
  -o OUTPUT, --output OUTPUT
                        Output file (default: output.csv)
```
This command will create an output.csv file in the current directory.

# Input and Output
The script requires the tracking file that was generated by StringTie. Usually, this file contains the extension .tracking.
The output csv file will have the transfrag ID, type of the transfrag, exon number, FPKM, FPKM-STD, TPM, TPM-STD, coverage, coverage STD, length, and length STD.

# Example of how to use this pipeline
All files used in this example can be found in the folder test.
```shell
python identical.py -t test/ALL_merged_track.tracking
                    -f test/ALL_merged_track.fasta
                    -l test/lncrna_ids.txt
                    -m [x,i,u,=]
                    -g g1:[1,2,3]/g2:[4,5,6]
```
- ALL_merged_track.tracking content:
```shell
TCONS_00000001	XLOC_000001	gene-ZemaCp004|rna-ZemaCp004	s	q1:STRG.34410|STRG.34410.1|2|4.936105|7.642557|39.146103|616	q2:STRG.33920|STRG.33920.1|2|2.702848|4.310839|22.575748|3225	q3:STRG.36866|STRG.36866.1|2|3.969663|5.404479|37.650860|626	q4:STRG.32721|STRG.32721.1|2|3.605668|5.286214|24.715002|616	q5:STRG.35714|STRG.35714.1|2|6.009348|7.791325|52.893368|615	q6:STRG.35954|STRG.35954.1|2|5.127873|5.959642|48.961872|584
TCONS_00000002	XLOC_000002	gene-ZemaCp005|rna-ZemaCp005	=	q1:STRG.34411|STRG.34411.1|1|1.285085|1.989693|10.191448|186	q2:STRG.33921|STRG.33921.1|1|0.652475|1.040649|5.449849|186	q3:STRG.36867|STRG.36867.1|1|0.494043|0.672612|4.685820|186	q4:STRG.32722|STRG.32722.1|1|0.470821|0.690263|3.227235|186	q5:STRG.35715|STRG.35715.1|1|1.023219|1.326638|9.006216|186	q6:STRG.35955|STRG.35955.1|1|0.786290|0.913831|7.507645|186
TCONS_00000003	XLOC_000003	gene-ZemaCp006|rna-ZemaCp006	=	q1:STRG.34412|STRG.34412.1|1|0.262803|0.406897|2.084178|111	-	-	-	q5:STRG.35716|STRG.35716.1|1|0.131300|0.170236|1.155687|111	q6:STRG.35956|STRG.35956.1|1|0.139318|0.161916|1.330236|111
TCONS_00000004	XLOC_000004	gene-ZemaCp008|rna-ZemaCp008	=	q1:STRG.34413|STRG.34413.1|1|0.351223|0.543797|2.785395|1422	q2:STRG.33922|STRG.33922.1|1|0.171253|0.273135|1.430401|1422	q3:STRG.36868|STRG.36868.1|1|0.198366|0.270064|1.881428|1422	q4:STRG.32723|STRG.32723.1|1|0.196361|0.287882|1.345955|1422	q5:STRG.35718|STRG.35718.1|1|0.290750|0.376968|2.559141|1422	q6:STRG.35957|STRG.35957.1|1|0.205374|0.238687|1.960947|1422
TCONS_00000005	XLOC_000005	gene-ZemaCp011|rna-ZemaCp011	=	q1:STRG.34416|STRG.34416.1|1|1.764536|2.732026|13.993769|105	q2:STRG.33924|STRG.33924.1|1|0.861831|1.374556|7.198512|105	q3:STRG.36870|STRG.36870.1|1|0.910000|1.238915|8.631027|105	q4:STRG.32725|STRG.32725.1|1|0.714448|1.047441|4.897172|105	q5:STRG.35720|STRG.35720.1|1|1.632656|2.116794|14.370388|105	q6:STRG.35960|STRG.35960.1|1|1.782282|2.071378|17.017559|105
```
- ALL_merged_track.fasta content:
```shell
>TCONS_00000001 loc:NC_001666.2|89-1150|- exons:89-1150 segs:1-1062 gene_name=psbA;oId=STRG.34408.1;cmp_ref=rna-ZemaCp002;class_code==;tss_id=TSS48;num_samples=6
ATGACTGCAATTTTAGAGAGACGCGAAAGTACAAGCCTGTGGGGTCGCTTCTGCAACTGGATAACTAGCA
CCGAAAACCGTCTTTACATTGGATGGTTCGGTGTTTTGATGATCCCTACCTTATTGACCGCAACTTCCGT
ATTTATTATCGCCTTCATCGCTGCTCCTCCAGTAGATATTGATGGTATTCGTGAGCCTGTTTCTGGTTCT
TTACTTTATGGAAACAATATTATCTCTGGTGCTATTATTCCTACTTCTGCGGCGATAGGGTTGCATTTTT
ACCCAATTTGGGAAGCTGCATCTGTTGATGAATGGTTATACAATGGCGGTCCTTATGAGCTAATTGTTCT
ACACTTCTTACTTGGTGTAGCTTGTTATATGGGTCGTGAGTGGGAACTTAGTTTCCGTCTGGGTATGCGC
CCTTGGATTGCTGTTGCATATTCAGCTCCTGTTGCAGCTGCTACTGCTGTTTTCTTGATTTACCCTATTG
GTCAAGGAAGTTTCTCTGATGGTATGCCTTTAGGAATATCTGGTACTTTCAACTTTATGATTGTATTCCA
GGCAGAGCACAACATCCTTATGCATCCATTTCACATGTTAGGTGTAGCTGGTGTATTCGGCGGTTCCCTA
TTCAGTGCTATGCATGGTTCCCTTGTAACCTCTAGTTTGATCAGGGAAACCACTGAAAATGAGTCTGCTA
ATGAGGGTTACAAATTTGGTCAGGAAGAAGAGACTTATAATATTGTGGCTGCTCACGGTTATTTTGGTCG
ATTAATCTTCCAATATGCTAGTTTCAACAATTCTCGTTCTTTACACTTCTTCTTGGCTGCTTGGCCTGTA
GTAGGGATCTGGTTCACTGCTTTAGGTATTAGTACTATGGCATTCAACCTAAATGGTTTCAATTTCAACC
AATCTGTAGTTGATAGCCAAGGTCGCGTTATTAATACTTGGGCTGATATCATCAACCGTGCTAATCTTGG
TATGGAGGTAATGCACGAACGTAATGCTCACAACTTCCCTCTAGACCTAGCTGCTCTTGAAGTTCCTTAC
CTTAATGGATAA
```

- GFF compare transfag class codes:
<p align="center">
  <img src="https://ccb.jhu.edu/software/stringtie/gffcompare_codes.png" width="350" alt="From GFF compare tool">
</p>
https://ccb.jhu.edu/software/stringtie/gffcompare.shtml#transfrag-class-codes

- lncrna_ids.txt content:
```shell
TCONS_00000016
TCONS_00000031
TCONS_00000020
TCONS_00000021
TCONS_00000026
TCONS_00000029
TCONS_00000034
TCONS_00000035
```

- ```shell
g1:[1,2,3]/g2:[4,5,6]
```
The g1 and g2 refer to the number of experiments you have performed using StringTie2. In this case, there are two experiments: control and treatment. Each experiment has 3 replicates, numbered 1, 2, and 3.

# Bugs
Do you have an issue? Please file it on :https://github.com/rogergsilva/structurally-identical-lncrnas/issues

# Citation
