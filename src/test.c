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
    char **b;
} B;

void test(bool t);

int main () {

    char *seq = "ATTAAAGGTTTATACCTTCCCAGAG";
    int length = strlen(seq);
    printf("%d\n", length);

    B *m = malloc(sizeof(B));
    m->b = (char **)bw(seq, length);
    printBWM(m->b, length);

    free(m->b);
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