#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <ctype.h>
#include "fasta.h"

/* takes a ref-sequence, builds fmIndex, and write to output */
int fmIndex();

/* takes the fmIndex of a ref-sequence, a file of reads, and aligns them */
int align();

/* runs align and fmIndex */
/* ./fmmap ref_seq indexOutput reads alignOutput */
int main(int argc, char **argv) {

    /* invalid number of arguments */
    if (argc != 5) {
        printf("invalid number of arguments.");
        return 0;
    }

    /* string path to each file */
    char *seq = argv[1];
    char *indexOut = argv[2];
    char *reads = argv[3];
    char *alignOut = argv[4];

    /* run fmIndex */
    fmIndex();

   return 0;
}

int fmIndex() {

    FASTAFILE *ffp;
    char *seq;
    char *name;
    int L;

    ffp = OpenFASTA("2019-nCov.fa");
    while (ReadFASTA(ffp, &seq, &name, &L)) {
        printf(">%s\n", name);
        printf("%s\n", seq);

        free(seq);
        free(name);
    }
    CloseFASTA(ffp);


    return 0;
}

int align() {
    return 0;
}