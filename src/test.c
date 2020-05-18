#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <ctype.h>
#include "fasta.h"
#include "auxiliary.h"

typedef struct sa {
    int *suffixArray;
} suffArr;

typedef struct B {
    char *name;
    char *seq;
    int length;
    int *suffixArray;
    char **b;
} B;

void test(bool t);

int main () {
/*
    char *seq = "ATTAAAGGTTTATACCTTCCCAGGTAACAAACCAACCAACTTTCGATCTCTTGTAGATCTGTTCTCTAAA";
    int length = strlen(seq);
    printf("%d\n", length);

    B *m = malloc(sizeof(B));

    m->b = (char **)bw(seq, length);
    printBWM(m->b, length);
    free(m);
*/

    /* lets try w fasta parser  */
    /* parse .fa file contains our reference sequence and build fm-index */
    char *reference = "./smallFasta.fa";
    FASTAFILE *ffp;
    char *seq;
    char *name;
    int length;

    ffp = OpenFASTA(reference);
    B *m = malloc(sizeof(B));
    while (ReadFASTA(ffp, &seq, &name, &length)) {

        /* name */
        m->name = malloc(strlen(name) + 1);
        strcpy(m->name, name);

        /* seq */
        m->seq = malloc(strlen(seq) + 1);
        strcpy(m->seq, seq);
        
        /* length of seq */
        m->length = length;

        /* suffix array */
        m->suffixArray = (int *)buildSuffixArray(seq, length);
        printSA(m->suffixArray, length);

        /* burrows-wheeler matrix (bwm) */
        m->b = (char **)bw(seq, length);
        printBWM(m->b, length);

        free(seq);
        free(name);
    }

    CloseFASTA(ffp);

    free(m->b);
    free(m->suffixArray);
    free(m->seq);
    free(m->name);
    free(m);


    return 0;
}

void test(bool t) {
    if (!t) {
        printf("false");
    } else {
        printf("true");
    }
}