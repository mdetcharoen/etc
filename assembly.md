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


# 1. Obtain fastq files
Please download fastq files from NCBI. We will try to assemble a genome of E. coli from Illumina reads.
https://trace.ncbi.nlm.nih.gov/Traces/?view=run_browser&acc=SRR27958868&display=download


## 2. Quality control using fastp
fastp is a all-in-one program to process fastq files.
install fasp via conda or mamba using

``
conda install fastp
``

or 

``
mamba install fastp
``


fastp -i SRR27958868.fastq.gz -o out_SRR27958868.fastq.gz





