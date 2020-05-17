#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <ctype.h>
#include <stdbool.h>

#include "fasta.h"
#include "fm.h"
#include "auxiliary.h"

/* takes a ref-sequence, builds fmIndex, and writes to output */
int fmIndex(char *reference, char *output);

/* takes the fmIndex of a ref-sequence, a .fa of reads, and aligns them */
int align();

/* runs align and fmIndex */
/* ./fmmap ref_seq indexOutput reads alignOutput */
int main(int argc, char **argv) {

    /* invalid number of arguments
    if (argc != 5) {
        printf("invalid number of arguments.");
        return 0;
    }
    */

    /* string path to each file */
    char *ref_seq = argv[1];
    char *indexOut = argv[2];
    // char *reads = argv[3];
    // char *alignOut = argv[4];

    /* fmIndex */
    fmIndex(ref_seq, indexOut);

   return 0;
}

/* reference is a path to a .fa containing a genetic-sequence */
int fmIndex(char *reference, char *output) {

    /* parse .fa file contains our reference sequence and build fm-index */
    FASTAFILE *ffp;
    char *seq;
    char *name;
    int length;

    ffp = OpenFASTA(reference);
    FM *fm = malloc(sizeof(FM));
    while (ReadFASTA(ffp, &seq, &name, &length)) {

        /* name */
        fm->name = malloc(strlen(name) + 1);
        strcpy(fm->name, name);

        /* seq */
        fm->seq = malloc(strlen(seq) + 1);
        strcpy(fm->seq, seq);
        
        /* length of seq */
        fm->length = length;

        /* suffix array */
        fm->suffixArray = (int*)(buildSuffixArray(seq, length, false));

        /* burrows-wheeler matrix (bwm) */


        /* occTable */

        free(seq);
        free(name);
    }

    CloseFASTA(ffp);

    /* free our fm-index */
    free(fm->seq);
    free(fm->name);
    free(fm);

    return 0;
}

int align() {
    return 0;
}