## <ins>COVID-19 Genome Variant, Burrows-Wheeler Transform, and FM-Index </ins>

## 

![10k Reads](./IGV-10K.png)

A read-alignment tool that utilizes the seed-and-extend paradigm, FM-Index, suffix arrays, burrows-wheeler transform, and core dynamic programming principles to provide for efficient seed-searching to discover variants in the COVID-19 genome.

## 

### (1) <ins>Clone</ins>
```
git clone https://github.com/RyanDBurke/coronavirus-genome-variant.git
```

##

### (2) <ins>Execution</ins>

##### first, compile
```
$ cd src/
$ make
```

##### then, execute one of the valid commands below
```
$ ./fmmap covid 1K
$ ./fmmap covid 10K
$ ./fmmap covid 1M
$ ./fmmap covid default
$ ./fmmap <reference-sequence>.fa <output file>.txt <reads>.fa.gz <output file>.sam
```

##

### (3) <ins>IGV Visualization</ins>
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
&nbsp;&nbsp;&nbsp;&nbsp;(3a) Load genome: <em>Genomes -> Load Genome From File</em> and select ```2019-nCoV.fa```
&nbsp;&nbsp;&nbsp;&nbsp;(3b) Load alignment: <em>File -> Load From File</em> and select ```mapping_sorted.bam```
&nbsp;&nbsp;&nbsp;&nbsp;optional: details on how to navigate IGV [here](https://software.broadinstitute.org/software/igv/AlignmentData)


### <ins>Structure</ins>
    coronavirus-genome-variant
    ├── README                   
    └── src
        ├── 2019-nCov.fa                # nCov-19 Genome
        ├── Default                     # default input/output
        |   ├── reads-small.fa          
        |   └── ref-small.fa
        |
        ├── FM-output                   # standard output file for FM-Index          
        |    └── FMindex.txt
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
        ├── fmmap.c                     # builds our FM-Index and Aligner
        ├── fmmap.h                     # header
        └── kseq.h                      # FASTA parser
##

### <ins>Auxiliary Links</ins>
* [FASTA parser](https://github.com/lh3/readfq)
* [nCov-19 Genome](https://www.ncbi.nlm.nih.gov/nuccore/NC_045512.2?report=fasta)

##

### <ins> Future Goals </ins>
* [Compressed suffix array](https://www.cs.cmu.edu/~dga/csa.pdf), for more efficient look-ups
* Fast-rank calculations on burrows-wheeler transform

## 
Written for my Bioinformatic Algorithms, Databases, and Tools course<br />
Shoutout Dr. Rob Patro for the cool project!
