#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <ctype.h>
#include "fasta.h"
#include "auxiliary.h"

typedef struct sa {
    int *suffixArray;
} suffArr;

void test();

int main () {

    char *seq = "banana";

    suffArr *s = malloc(sizeof(suffArr));

    s->suffixArray = (int*)(buildSuffixArray(seq, 6));

    // int *sa = buildSuffixArray(seq, 6);

    printSA(s->suffixArray, 6);



    /*
    for (int i = 0; i < 6; i++) {
        printf("%d\n", s->suffixArray[i]);
    }
    */

   free(s->suffixArray);
    free(s);

    return 0;
}