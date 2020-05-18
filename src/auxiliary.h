/* auxiliary.h:  contains all aux-functions needed to build index and align reads */

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <ctype.h>
#include <stdbool.h>

int *buildSuffixArray(char *seq, int length);
char **bw(char *seq, int length);
void printSA(int *sa, int length);
void printBWM(char **b, int length);

/* suffix array struct */
typedef struct suffixArray {
    int offset;
    char *suffix;
} SA;

/* rotation struct for bwm */
typedef struct R {
    int offset;
    char *rotation;
} R;

/* comparison sort function for suffix arrays */
int cmpSA(const void *a, const void *b) {
    const SA *da = (const SA *) a;
    const SA *db = (const SA *) b;

    return strcmp(da->suffix, db->suffix);
}

/* comparison sort function for bwm */
int cmpBMW(const void *a, const void *b) {
    const R *da = (const R *) a;
    const R *db = (const R *) b;

    return strcmp(da->rotation, db->rotation);
}

/* builds suffix array */
int *buildSuffixArray(char *seq, int length) {

    /* store suffixes and their offset */
    SA suffixes[length];

    /* starting adding SA objects to array */
    for (int i = 0; i < length; i++) {
        suffixes[i].suffix = malloc(strlen(seq + i) + 1);
        strcpy(suffixes[i].suffix, (seq + i));
        suffixes[i].offset = i;
    }

    /* sort suffixes */
    qsort(suffixes, length, sizeof(SA), &cmpSA);

    /* build suffix array */
    int *sa = malloc(sizeof(int) * length);
    for (int i = 0; i < length; i++) {
        sa[i] = suffixes[i].offset;
        free(suffixes[i].suffix);
    }

    return sa;
}

/* build burrows-wheeler matrix */
char **bw(char *seq, int length) {

    /* store rotations and their offset */
    R bwMatrix[length];

    /* starting adding rotations */
    for (int i = 0; i < length; i++) {

        /* offset */
        bwMatrix[i].offset = i;

        /* build rotation */
        char *prefix = malloc(i + 1);
        bwMatrix[i].rotation = malloc(length);
        strcpy(bwMatrix[i].rotation, (seq + i));
        strncpy(prefix, seq, i);
        strcat(bwMatrix[i].rotation, prefix);

        /* free used memory */
        free(prefix);
    }

    /* sort rotations */
    qsort(bwMatrix, length, sizeof(R), &cmpBMW);

    /* stores the actual matrix, with just the rotations */
    char **matrix = malloc(sizeof(char *) * length);

    /* fill matrix */
    for (int i = 0; i < length; i++) {
        matrix[i] = malloc(sizeof(char) * length);
        strcpy(matrix[i], bwMatrix[i].rotation);
        free(bwMatrix[i].rotation);
    }
    
    return matrix;
}

/* prints suffix array */
void printSA(int *sa, int length) {
    for (int i = 0; i < length; i++) {
        printf("%d\n", sa[i]);
    }
}

/* prints bwm */
void printBWM(char **b, int length) {
    for (int i = 0; i < length; i++) {
        printf("%s\n", b[i]);
    }
}

