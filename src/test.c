#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <ctype.h>
#include "fasta.h"

void test();

int main () {
    
    int *num = malloc(sizeof(int*) * 4);
    test(num);
    int i;
    for (i = 0; i < 10; i++) {
        printf("%d\n", num[i]);
    }

    free(num);

    return 0;
}

void test(int *a) {

    int i;
    for (i = 0; i < 10; i++) {
        a[i] = i;
    }
    return a;
}