#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <ctype.h>
#include <stdbool.h>

#include "fmmap.h"

/*

    DONT DELETE!!!!
    this was great for debugging bwm and suffix.
    you found prefixes need to be null terminated
    and a an array of struct with char* dont need
    to be free'd if you used malloc (idky)

*/

typedef struct t {
    char *suffix;
} T;

char **test();
void print(char**);
int cmpT(const void *a, const void *b);

int main() {

    char *seq = "ATTAAAGGTTTATTAAAGGTTTATTAAAGGTTTATTAAAGGTTTATTAAAGGTTTATTAAAGGTTTATTAAAGGTTTATTAAAGGTTTATTAAAGGTTTATTAAAGGTTTATTAAAGGTTTATTAAAGGTTTATTAAAGGTTTATTAAAGGTTTATTAAAGGTTT";
    int length = strlen(seq);

    printf("length: %d\n\n", length);

    T array[length];

    for (int i = 0; i < length; i++) {
        char *prefix = malloc(i + 1);

        printf("i: %d\n", i);

        char *suff = (seq + i);
        array[i].suffix = malloc(length + 1);
        strcpy(array[i].suffix, suff);

        printf("suffix: %s\n", array[i].suffix);
        printf("memory allocated for result: %d\n", (length + 1));

        strncpy(prefix, seq, i);

        printf("prefix: %s\n", prefix);
        printf("memory allocated for prefix: %d\n", (i + 1));
        strcat(array[i].suffix, prefix);

        printf("Result: %s\n\n\n", array[i].suffix);

        free(prefix);
    }

    /* sort suffixes */
    qsort(array, length, sizeof(T), &cmpT);

    /* debugging purposes */
    printf("Matrix sorted:\n");
    for (int i = 0; i < length; i++) {
        printf("%s\n", array[i].suffix);
    }

    int memory = sizeof(char *) * (length + 1);
    char **matrix = malloc(sizeof(char *) * (length + 1));
    printf("\n\nmemory allocated for matrix: %d\n", memory);

    /* fill matrix */
    printf("matrix actual:\n");
    for (int i = 0; i < length; i++) {
        matrix[i] = malloc(sizeof(char) * length);
        strcpy(matrix[i], array[i].suffix);
        printf("%s\n", matrix[i]);
        // free(array[i].suffix);
    }

    free(matrix);


    return 0;
}

/* comparison sort function for bwm */
int cmpT(const void *a, const void *b) {
    const T *da = (const T *) a;
    const T *db = (const T *) b;

    return strcmp(da->suffix, db->suffix);
}

void print(char **s) {
    for (int i = 0; i< 10; i++) {
        printf("%s\n", s[i]);
    }
}

