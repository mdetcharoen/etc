# Bacterial genome assembly

## Why assemble a genome
To determine the complete genome sequence of an organism/virus.

## Steps in genome assembly
1. Obtain sequence read file(s).
2. Quality control.
3. Raw data cleanup/quality trimming.
4. Assembly using appropiate assembler into conigs/scaffolds. Short-read assembler, long-read assembler, or long+short-read assembler.
5. Examine the output of the assembly and assess assembly quality.


Reads can be saved in a Fasta file as text or in a FastQ file along with their qualities. They can also be saved as alignments to references in different formats like SAM or its binary compressed version, BAM. Except for the binary BAM format, all other file formats can be compressed and are often stored in a compressed form (with .gz extension for gzipped files).

-----

### Install SRA toolkit
```
wget --output-document sratoolkit.tar.gz https://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/3.0.10/sratoolkit.3.0.10-ubuntu64.tar.gz
```

```
tar -vxzf sratoolkit.tar.gz
```

```
export PATH=$PATH:$PWD/sratoolkit.3.0.10-ubuntu64/bin
```

```
which fastq-dump
```

This should return something like: ``/home/matt/miniconda3/bin/fastq-dump``

<br/>

### add bioconda channel to your conda

copy these commands into your terminal and hit enter

```
conda config --add channels defaults
conda config --add channels bioconda
conda config --add channels conda-forge
conda config --set channel_priority strict
```
-----
<br/>

## 1. Obtain fastq files
Please download fastq files from NCBI. We will try to assemble a genome of E. coli from Nanopore reads.

https://trace.ncbi.nlm.nih.gov/Traces/?view=run_browser&page_size=10&acc=SRR27502765&display=download


Count how many sequences in the fastq file.

<br/>


for unzipped fastq file:

```
echo $(cat SRR27502765.fastq|wc -l)/4|bc
```


for .gz file:

```
echo $(zcat SRR27502765.fastq.gz|wc -l)/4|bc
```



<br/>

## 2. Quality control 

### using FastQC

FastQC (https://www.bioinformatics.babraham.ac.uk/projects/fastqc/) is a basic program to check your sequencing quality.

<br/>

### using fastp
fastp (https://github.com/OpenGene/fastp) is a all-in-one program to process fastq files.
install fasp via conda or mamba using

```
conda install fastp
```

or 

```
mamba install fastp
```

check the instalation with a command

```
fastp
```

version?

```
fastp -v
```

list all command using 

```
fastp -h
```

Have a look at ``-M`` function. 

``-M, --cut_mean_quality`` the mean quality requirement option shared by cut_front, cut_tail or cut_sliding. Range: 1~36 default: 20 (Q20) (int [=20])



``-l``

``--length_required``
reads shorter than length_required will be discarded, default is 15. (int [=15])

``--length_limit``
reads longer than length_limit will be discarded, default 0 means no limitation. (int [=0])

<br/>
<br/>


Now check the reads we just downloaded. 

```
fastp -i SRR27502765.fastq.gz -o outSRR27502765.fastq.gz
```

Fastp will produce 2 output files, fastp.json and fastp.html. 

How many reads before and after the filtering.

If we want to keep reads longer than 250 bp with quality >= Q25, what fastp command should you use?


Although we will not use the original file anymore, never delete it. 

<br/>

## 3. Assembly using Flye

Flye (https://github.com/fenderglass/Flye) is a popular long-read genome assembler. Please download it via conda/mamba

```
mamba install flye
```

check flye version

```
flye -v
```

list all commands

```
flye -h
```


Here we have the filtered Nanopore data. We know that E. coli genome should be around 4-6 Mb.

```
flye --nano-raw outSRR27502765.fastq.gz --out-dir outflye --thread 12 -g 6m -m 1000
```


<br/>

## 4. Check the assembly

### QUAST

Get basic statistics from QUAST (https://quast.sourceforge.net/quast.html)

```
conda install quast
```

```
quast -h
```

<br/>

now check the assembly:


```
quast assembly.fasta -t 12
```

