/* auxiliary.h:  contains all aux-functions needed to build index and align reads */

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <ctype.h>
#include <stdbool.h>

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

int *buildSuffixArray(char *seq, int length);
char **bw(char *seq, int length);
void printSA(int *sa, int length);
void printBWM(char **b, int length);

