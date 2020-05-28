## <ins>COVID-19 Genome Variant, Burrows-Wheeler Transform, and FM-Index </ins>

## 

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
```
```
$ ./fmmap covid 10K
```
```
$ ./fmmap covid 1M
```
```
$ ./fmmap covid default
```
```
$ ./fmmap <reference-sequence>.fa <output file>.txt <reads>.fa.gz <output file>.sam
```

##

### (3) <ins>IGV Visualization</ins>
##### install samtools
```
$ sudo apt-get update -y
$ sudo apt-get install -y samtools
```
##### install [IGV](https://software.broadinstitute.org/software/igv/download)
##### go-to mapping
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


### <ins>Structure</ins>
    BWT
    ├── README                   
    └── src
        ├── Makefile                    # build executables
        ├── fmmap.c                     # builds our FM-Index and Aligner
        ├── fmmap.h                     # header
        ├── ref_CoV19.fa                # nCov-19 Genome
        ├── reads.fa.gz                 # simulated coronavirus variant reads
        └── index_out                   # output of FM-Index built)  
##

### <ins>Auxiliary Links</ins>
* [FASTA parser](https://github.com/lh3/readfq)
* [nCov-19 Genome](https://www.ncbi.nlm.nih.gov/nuccore/NC_045512.2?report=fasta)

##

### <ins> Future Goals </ins>
* [Compressed suffix array](https://www.cs.cmu.edu/~dga/csa.pdf), for more efficient look-ups
* Fast-rank calculations on burrows-wheeler transform

## 
Written for my Bioinformatic Algorithms, Databases, and Tools course
