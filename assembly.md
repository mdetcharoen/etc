# Introduction to Genome Assembly

Genome assembly is a crucial step in genomics to reconstruct the complete genome sequence of an organism from fragmented DNA sequences. This process allows researchers to understand the genetic blueprint of an organism, identify genes, and explore its functional biology.

There are two primary approaches to genome assembly:

1. **Reference-Based Assembly**: This method aligns sequencing reads to a known reference genome. It is computationally efficient and suitable when a closely related reference genome is available. However, it may miss novel sequences not present in the reference genome.

2. ***de novo* Assembly**: This method assembles reads without relying on a reference genome. It is more computationally intensive but essential for studying organisms without a close reference or to identify unique sequences.

---

### Why Assemble a Genome?

- **Comprehensive Understanding**: To determine the complete genetic composition of an organism or virus.
- **Novel Discovery**: Identification of novel genes, operons, and regulatory elements.
- **Applications**: Used in microbial studies, evolutionary biology, medicine, agriculture, and environmental sciences.

---

### Steps in Genome Assembly

1. **Obtain Sequence Read Files**: Obtain sequencing data in FASTQ format from public repositories like NCBI.

2. **Quality Control**:
   - Check read quality using tools like **FastQC**.
   - Perform trimming and filtering with tools like **fastp** to improve the data quality.

3. **Raw Data Cleanup/Quality Trimming**:
   - Remove low-quality bases, adapters, and artifacts.
   - Use specified quality thresholds (e.g., Q20, Q25).

4. **Assembly**:
   - Choose an appropriate assembler:
     - Short-read assemblers (e.g., SPAdes for Illumina reads).
     - Long-read assemblers (e.g., Flye for Nanopore or PacBio reads).
     - Hybrid assemblers for combined datasets.
   - Assemble reads into contigs and scaffolds.

5. **Assess Assembly Quality**:
   - Use tools like **QUAST** to evaluate contig length, N50 values, and assembly accuracy.
   - Check for completeness with tools like **BUSCO**.

6. **Annotation**:
   - Annotate genes using tools like **Prokka**.

7. **Visualization and Mapping**:
   - Visualize the assembly graph with tools like **Bandage**.
   - Map reads to a reference using alignment tools like **minimap2**.

---

## **Practical Steps**

#### Step 1: Install Required Tools
Add **bioconda** to your Conda environment:
```bash
conda config --add channels defaults
conda config --add channels bioconda
conda config --add channels conda-forge
conda config --set channel_priority strict
```

---

### **Genome Assembly Practical**

#### **Step 2: Download Sequencing Data**
Download an *E. coli* dataset from NCBI:
[Link to Data](https://trace.ncbi.nlm.nih.gov/Traces/?view=run_browser&page_size=10&acc=SRR27502765&display=download).

Use the command below to download and check how many sequences are in the FASTQ file.

- If the file is **unzipped**:
  ```bash
  echo $(cat SRR27502765.fastq | wc -l) / 4 | bc
  ```
- If the file is **gzipped**:
  ```bash
  echo $(zcat SRR27502765.fastq.gz | wc -l) / 4 | bc
  ```

---

#### **Step 3: Quality Control**

**Install [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)**:
```bash
conda install fastqc
```

Run FastQC:
```bash
fastqc SRR27502765.fastq.gz
```

Generate a summary of sequencing quality.

**Install [fastp](https://github.com/OpenGene/fastp)**:
```bash
conda install fastp
```

Run fastp for filtering and trimming:
```bash
fastp -i SRR27502765.fastq.gz -o filtered.fastq.gz
```


Fastp will produce 3 output files, one filtered fastq file and two reports, fastp.json and fastp.html. 


![image](https://github.com/mdetcharoen/etc/assets/70691598/3178f2b4-bf63-4c78-9847-0c801d8876c3)


<br/>

How many reads before and after the filtering.


<details> 
  <summary>Fastp commands</summary>
  
 - list all command using 
  ```
  fastp -h
  ```

 - Have a look at ``-M`` function. 

  ``-M, --cut_mean_quality`` the mean quality requirement option shared by cut_front, cut_tail or cut_sliding. Range: 1~36 default: 20 (Q20) (int [=20])

  
 - reads shorter than length_required will be discarded, default is 15. (int [=15])
  ``-l,--length_required``
 
 - reads longer than length_limit will be discarded, default 0 means no limitation. (int [=0])
  ``--length_limit``

</details>

</br>
If we want to keep reads longer than 250 bp with quality >= Q25, what fastp command should you use?

<details>
  <summary>answer</summary>
  To keep reads longer than 250 bp and quality ≥ Q25:
  
```bash
fastp -i SRR27502765.fastq.gz -o filtered.fastq.gz -l 250 -M 25
```

</details>


Although we will not use the original file anymore, never delete it. 

---

#### **Step 4: Assembly**

**Install Flye**:
```bash
mamba install flye
```

Here we have the filtered Nanopore data. We know that E. coli genome should be around 4-6 Mb. We set minimum overlap between reads to 1000 base pairs.

Run Flye for genome assembly:
```bash
flye --nano-raw filtered.fastq.gz --out-dir flye_output --threads 12 -g 6m
```


<details> 
  <summary>Flye output files </summary>
  
  1. assembly.fasta (https://github.com/mdetcharoen/etc/blob/main/assembly.fasta)

  2. assembly_graph.gfa (https://github.com/mdetcharoen/etc/blob/main/assembly_graph.gfa)

  3. assembly_graph.gv (https://github.com/mdetcharoen/etc/blob/main/assembly_graph.gv)

  4. assembly_info.txt (https://github.com/mdetcharoen/etc/blob/main/assembly_info.txt)

  5. flye.log (https://github.com/mdetcharoen/etc/blob/main/flye.log)
     
</details>

<br/>



---

#### **Step 5: Evaluate Assembly**

**assembly graph**

We will use [Bandage](https://rrwick.github.io/Bandage/) to see the assembly graph.

Open Bandage > File > Load graph > select .gfa file (assembly_graph.gfa) > Draw graph

![image](https://github.com/mdetcharoen/etc/assets/70691598/12e4a5a6-ad48-46b0-904c-32b279b9b2b4)


**QUAST**:
```bash
conda install quast
```

Run QUAST to evaluate the assembly:
```bash
quast flye_output/assembly.fasta -t 12
```

**BUSCO**:
Install BUSCO:
```bash
conda install busco
```

Check gene completeness:
```bash
busco -i flye_output/assembly.fasta -m genome --auto-lineage-prok
```

---

#### **Step 6: Gene Annotation**

**Prokka**:
```bash
conda install prokka
```

Annotate the genome:
```bash
prokka --outdir annotated_genome --prefix ecoli flye_output/assembly.fasta
```

---
#### **Step 7: Mapping raw reads to the assembled genome**

**Install Minimap2**:
```bash
conda install minimap2
```

Map reads to the reference genome:
```bash
minimap2 -ax map-ont assembly.fasta SRR27502765.fastq.gz -o aligned.sam
```

**Convert SAM to BAM**:
Install Samtools:
```bash
conda install samtools
```

Convert and sort:
```bash
samtools sort aligned.sam > aligned.sorted.bam
```

Index the BAM file:
```bash
samtools index aligned.sorted.bam
```

---
#### **Step 8: Visualize the Results**

Install **IGV** (Integrative Genomics Viewer):
- [Download IGV](https://igv.org).

Load the `.bam` file in IGV to visualize the mapping results.

---

### **Explanation of Key Commands and Formats**

- **FASTQ**: Contains raw reads and their quality scores.
- **FASTA**: Simple sequence data without quality scores.
- **SAM/BAM**: Formats for storing alignments (SAM is text-based, BAM is binary).
- **N50**: A metric to assess assembly contiguity.
- **gfa**: A file format for genome assembly graphs.

---
















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
