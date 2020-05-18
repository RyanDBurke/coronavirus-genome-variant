/* fm.h */

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <ctype.h>

#define MAX_LEN 40000

typedef struct fm {
    char *name;
    char *seq;
    int length;
    int *suffixArray;
    char **bwm;
    // occTable
} FM;