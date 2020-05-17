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

    char *s = "ryan";
    char *a = malloc(4 + 1);
    char *b = malloc(4 + 1);

    strcpy(a, (s + 1));
    strncpy(b, s, 1);
    strcat(a, b);

    free(b);
    free(a);

    printf("%s\n", a);

    
    return 0;
}

void test(bool t) {
    if (!t) {
        printf("false");
    } else {
        printf("true");
    }
}