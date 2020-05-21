###### Unfinished 5/21/20.

## <ins>COVID-19 Genome Variant, Burrows-Wheeler Transform, and FM-Index </ins>

## 

A read-alignment tool that utilizes the seed-and-extend paradigm, FM-Index, suffix arrays, and the burrows-wheeler matrix
to provide for efficient seed-searching to discover variants in the COVID-19 genome.

## 

### <ins>Execution</ins>

#### Default Execution
```
$ cd src/
$ make
$ ./fmmap default
```
#### Example Execution
```
$ cd src/
$ make
$ ./fmmap ref-CoV19.fa index_out.txt reads.fa.gz alignments.sam
```
#### General Execution
```
$ cd src/
$ make
$ ./fmmap <reference-sequence>.fa <index output file> <reads>.fa <align output file>.sam
```

##

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
* [FASTA parser](https://github.com/eturro/mmseq/blob/master/src/fasta.c)
* [nCov-19 Genome](https://www.ncbi.nlm.nih.gov/nuccore/NC_045512.2?report=fasta)

##

### <ins> Future Goals </ins>
* [Compressed suffix array](https://www.cs.cmu.edu/~dga/csa.pdf), for more efficient look-ups

## 
Written for my Bioinformatic Algorithms, Databases, and Tools course
