# <ins>Burrows-Wheeler Transform and COVID-19 Genome Variant </ins>

## 

> Implements string macthing using seed-and-extend strategy 
> on exact matching using the FM-Index and semi-global read-aligning.

## 

### <ins>Structure</ins>
    BWT
    ├── README                   
    └── src
        ├── fmmap.c                     (builds our FM-Index and Aligner)
        ├── auxiliary.h                 (auxiliary functions)
        ├── fasta.h                     [FASTA parser](https://github.com/eturro/mmseq/blob/master/src/fasta.c)
        ├── ref_CoV19.fa                [nCov-19 Genome](https://www.ncbi.nlm.nih.gov/nuccore/NC_045512.2?report=fasta)
        └── index_out                   (output of FM-Index built)   
##

### <ins> Future Goals </ins>
* [Compressed suffix array](https://www.cs.cmu.edu/~dga/csa.pdf), for more efficient look-ups

