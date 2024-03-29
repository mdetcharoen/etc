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

Fastp will produce 3 output files, one filtered fastq file and two reports, fastp.json and fastp.html. 


![image](https://github.com/mdetcharoen/etc/assets/70691598/3178f2b4-bf63-4c78-9847-0c801d8876c3)


<br/>

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



<details> 
  <summary>flye output files </summary>
  
  1. assembly.fasta (https://github.com/mdetcharoen/etc/blob/main/assembly.fasta)

  2. assembly_graph.gfa (https://github.com/mdetcharoen/etc/blob/main/assembly_graph.gfa)

  3. assembly_graph.gv (https://github.com/mdetcharoen/etc/blob/main/assembly_graph.gv)

  4. assembly_info.txt (https://github.com/mdetcharoen/etc/blob/main/assembly_info.txt)

  5. flye.log (https://github.com/mdetcharoen/etc/blob/main/flye.log)
     
</details>

<br/>

## 4. Check the assembly

### assembly graph

We will use Bandage (https://rrwick.github.io/Bandage/) to see assembly graph.

Open Bandage > File > Load graph > select .gfa file (assembly_graph.gfa) > Draw graph

![image](https://github.com/mdetcharoen/etc/assets/70691598/12e4a5a6-ad48-46b0-904c-32b279b9b2b4)


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

### BUSCO

Check gene set completeness using BUSCO (https://busco.ezlab.org/busco_userguide.html#installation-with-conda)

```
busco -i assembly.fasta -m genome --auto-lineage-prok
```

## 5. Annotation

Prokka (https://github.com/tseemann/prokka)

```
prokka --outdir mygenome file.fasta
```




-------

# Map fastq reads to a reference genome

There are popular mapping tools such as Bowtie (https://bowtie-bio.sourceforge.net/manual.shtml), BWA (https://github.com/lh3/bwa), mummer (https://github.com/mummer4/mummer), and minimap2 (https://github.com/lh3/minimap2).

We will use minimap2 to map Nanopore reads to Drosophila melanogaster mitochondrial genome.

Download Nanopore reads SRR25018161 from NCBI SRA (https://trace.ncbi.nlm.nih.gov/Traces/?view=run_browser&acc=SRR25018161&display=download)

and the reference genome NC_024511.2 (https://www.ncbi.nlm.nih.gov/nuccore/NC_024511.2?report=fasta)

1. Quality check the Nanopore data.

2. map to the reference using miimap2

```
minimap2 -ax map-ont -I8g mtDNA_Dmel_sequence.fasta wmel_SRR25018161.fastq.gz -o aligned_mtDNA.sam
```

How many reads mapped?

what is `.sam` file?

3. convert and sort `.sam` to `.bam` file

```
samtools sort aligned_mtDNA.sam > aligned_mtDNA.sorted.bam
```

what is `.bam` file?

4. index the `.sorted.bam` file

```
samtools index aligned_mtDNA.sorted.bam
```

5. check the alignment process using IGV (https://igv.org/doc/desktop/#DownloadPage/)
