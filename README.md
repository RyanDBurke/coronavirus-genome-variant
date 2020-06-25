# Table of Contents

* [What is This?](#what)
* [What Does it Do?](#cool)
* [How Do I Execute?](#execute)
* [File Structure](#structure)
* [Auxiliary Links](#links)
* [Future Goals](#goals)
* [Concluding Thoughts](#thoughts)

## What is This? <a name="what"></a>

An efficient seed-search alignment tool to compare read-variants to the COVID-19 virus genome using Burrows-Wheeler transform, suffix arrays, the FM-Index, and core dynamic programming principles.

## What Does it Do? <a name="cool"></a>

<figure>
  <img src="./IGV-10K.png" alt="Integrative Genomics Viewer for 10k reads" name="figure1">
</figure>

Steps:
1. Take the nCov-19 genome sequence and build the auxiliary structures (Burrows-Wheeler transform and suffix arrays) for the FM-Index
2. Take the 100bp-length simulated reads and align them to nCov-19 using the seed-and-extend paradigm and dynamic programming
3. Build the .SAM file mapping
4. Upload the mapping to Integrative Genomics Viewer (IGV)
5. View coverage track for an alignment ([Figure 1](#figure1))

## How Do I Execute? <a name="execute"></a>

### (1) Clone
```
git clone https://github.com/RyanDBurke/coronavirus-genome-variant.git
```

### (2) Execution

##### first, install zlib
```
sudo apt-get install libz-dev
```

##### then, compile
```
$ cd coronavirus-genome-variant
$ cd src/
$ make
```

##### then, execute ONE of the valid commands below
```
$ ./fmmap covid 1K
$ ./fmmap covid 10K
$ ./fmmap covid 1M
$ ./fmmap covid default
$ ./fmmap <reference-sequence>.fa <output file>.txt <reads>.fa.gz <output file>.sam
```

### (3) IGV Visualization
##### install samtools
```
$ sudo apt-get update -y
$ sudo apt-get install -y samtools
```
##### install IGV [here](https://software.broadinstitute.org/software/igv/download), then go-to your mapping
```
$ cd Mappings/
```
##### format for IGV
```
$ samtools view -b -o mapping.bam mapping.sam
$ samtools sort -o mapping_sorted.bam mapping.bam
$ samtools index mapping_sorted.bam
```
##### open IGV
&nbsp;&nbsp;&nbsp;&nbsp;**(3a)** Load genome: <em>Genomes -> Load Genome From File</em> and select ```2019-nCoV.fa``` <br />
&nbsp;&nbsp;&nbsp;&nbsp;**(3b)** Load alignment: <em>File -> Load From File</em> and select ```mapping_sorted.bam``` <br />
&nbsp;&nbsp;&nbsp;&nbsp;**optional**: details on how to navigate IGV [here](https://software.broadinstitute.org/software/igv/AlignmentData)

## File Structure <a name="structure"></a>
    coronavirus-genome-variant
    ├── README                   
    └── src
        ├── 2019-nCov.fa                # nCov-19 Genome
        ├── Default                     # default input/output
        |   ├── reads-small.fa          
        |   └── ref-small.fa
        |
        ├── FM-output                   # standard output file for FM-Index          
        |   └── FMindex.txt
        |
        ├── Makefile                    # build executables
        ├── Mappings                    # standard output .sam file(s) for alignments
        |   ├── mappings.sam                      
        |   └── mapping1M.sam           # holds the mapping for the 1M reads (so you don't have to wait 3hrs)
        |    
        ├── Reads                       # .gz FASTA files of reads (1K, 10K, and 1M)
        |   ├── reads_10K.fa.gz
        |   ├── reads_1K.fa.gz
        |   └── reads_1M.fa.gz
        |
        ├── fmmap.c                     # .c -> builds our FM-Index and Aligner
        ├── fmmap.h                     # .h -> struct/function-delcarations
        └── kseq.h                      # FASTA parser

## Auxiliary Links <a name="links"></a>
* [FASTA parser](https://github.com/lh3/readfq)
* [nCov-19 Genome](https://www.ncbi.nlm.nih.gov/nuccore/NC_045512.2?report=fasta)
* Execution guide for other OS distributions


## Future Goals <a name="goals"></a>
* [Compressed suffix array](https://www.cs.cmu.edu/~dga/csa.pdf), for more efficient look-ups
* Fast-rank calculations on burrows-wheeler transform

##### Concluding Thoughts<a name="thoughts"></a>
```
This project was originally given to us mid-semester, right when the chaos following the covid pandemic occurred.
I originally wrote it in python, and never finished it -- the scope of the assignment was beyond me and I was 
struggling to wrap my head around a lot of the intense algorithms.

After the end of the semester I decided to revisit the project using a compiled language (C) for faster 
runtimes and also for some well-need practice with C. I'm proud to say I finished the project and learned 
a TON about dynamic-programming, data-compression for efficient look-ups, file-structure manangement, 
memory/pointers, and plenty more.

All thanks goes out to my professor, Dr. Rob Patro, he was extremely helpful and I genuinely feel 
confident in my abilities as a programmer.
```
## 
<em>Written for my Bioinformatics course</em><br />
