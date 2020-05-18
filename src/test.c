#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <ctype.h>
#include "fasta.h"
#include "auxiliary.h"

typedef struct sa {
    int *suffixArray;
} suffArr;

void test(bool t);

int main () {

    char *seq = "banana$";
    char **m = (char **)bw(seq, 7);

    for (int i = 0; i < 7; i++) {
        printf("%s\n", m[i]);
    }

    free(m);


    /* allocate memory for number of strings
    char **b = malloc(3);

    for (int i = 0; i < 3; i++) {
        b[i] = malloc (50 * sizeof(char));
        strcpy(b[i], "hey");
        printf("%s\n", b[i]);
    }

    */
    
    return 0;
}

void test(bool t) {
    if (!t) {
        printf("false");
    } else {
        printf("true");
    }
}