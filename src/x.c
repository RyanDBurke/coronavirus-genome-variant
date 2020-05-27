#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <ctype.h>
#include <stdbool.h>
#include <math.h>

#include "fmmap.h"

typedef struct E {
    int score;
    int direction;
} Entry;

int test(int x);
void change(int *a);
void susbtring(char *res, char* string, int start, int end);
int intArray(int *a, int x, int y);
int maxAlign(int a, int b, int c);
int score(char a, char b, int gap);
int editDistance(int OPT[MAXROW][MAXCOL], char *x, char *y, int n, int m, int gap);
int minAlign(int a, int b, int c);
void printMatrix(int matrix[MAXROW][MAXCOL], int n, int m, char *x, char *y);

int main(int argc, char **argv) {

    char *x = "AAGGTATGAATCAA";
    char *y = "AAGGTATGAATC";
    int n = strlen(x) + 1;
    int m = strlen(y) + 1;
    int matrix[MAXROW][MAXCOL];
    int gap = 3;
    
    int edit = editDistance(matrix, x, y, n, m, gap);

    printf("edit distance: %d\n", edit);

    printMatrix(matrix, n, m, x, y);

    char *cigar = buildCigar(matrix, n, m, gap, x, y);

    printf("CIGAR: %s\n", cigar);
    
    free(cigar);

    return 0;

}

/* returns CIGAR string  */
char *buildCigar(int OPT[MAXROW][MAXCOL], int n, int m, int gap, char *x, char *y) {

    /* n and m are lengths, so we want < n and < m */
    n = n - 1;
    m = m - 1;

    /* result array we'll return */
    int traceIndex = 0;


    char *cigar = malloc(n * m);
    while (true) {

        /* we've reached matrix[0][0] */
        /* for fitting alignment it ends when we reach row = 0 */
        if (n == 0 || m == 0) {
            break;
        }

        int left = gap + OPT[n - 1][m];
        int diagonal = score(x[n - 1], y[m - 1], gap) + OPT[n - 1][m - 1];
        int down = gap + OPT[n][m - 1];

        int max = maxAlign(left, diagonal, down);
        if (max == diagonal) {
            n = n - 1;
            m = m - 1;
            cigar[traceIndex] = 'M';
            traceIndex++;
        } else if (max == left) {
            n = n - 1;
            cigar[traceIndex] = 'D';
            traceIndex++;
        } else {
            m = m - 1;
            cigar[traceIndex] = 'I';
            traceIndex++;
        }

    }


    return cigar;
}

int test(int x) {
    
    int a = 11;
    if (a > x) {
        x = a;
    }

    return x;
}

int min(int a, int b) {
    return (a > b) ? b : a;
}
int max(int a, int b) {
    return (a > b) ? a : b;
}

/* max function needed to fill OPT-matrix */
int maxAlign(int a, int b, int c) {
    return max(a, max(b, c));
}

/* min function for OPT-matrix */
int minAlign(int a, int b, int c) {
    return min(a, min(b, c));
}

/* return score between two characters */
int score(char a, char b, int gap) {
    if (a == b) {
        return 0;
    } else if (a != b) {
        return -1;
    } else if (a == '-'  || b == '-') {
        return gap;
    } else {
        return 0;
    }
}

void printMatrix(int matrix[MAXROW][MAXCOL], int n, int m, char *x, char *y) {

    for (int i = 1; i < m; i++) {
        if (i < 2) {
            printf("\t");
        } else {
            printf("\t%c", y[i - 2]);
        }
    }

    printf("\n\t");

    for (int i = 0; i <= n; i++) {
        if (i >= 1) {
            printf("%c\t", x[i]);
        }
        
        for (int j = 0; j <= m; j++) {
            printf("%d\t", matrix[i][j]);
        }
        printf("\n");
    }
}

/* builds OPT-matrix, and returns OPT[n][m] -- which is the edit distance between x and y */
int editDistance(int OPT[MAXROW][MAXCOL], char *x, char *y, int n, int m, int gap) {

    /* add our initial gap penalties to the first column of each row */
    for (int i = 0; i < n; i++) {
        OPT[i][0] = 0;
    }

    /* add our initial gap penalties to the first row of each column */
    for (int j = 0; j < m; j++) {
        OPT[0][j] = j * gap;
    }

    for (int i = 1; i <= n; i++) {
        for (int j = 1; j <= m; j++) {
            OPT[i][j] = maxAlign(
                    (score(x[i - 1], y[j - 1], gap) + OPT[i - 1][j - 1]),
                    (gap + OPT[i - 1][j]),
                    (gap + OPT[i][j - 1])
                ); /* calculate score */
        }
    }

    return OPT[n][m];
}

/* lowercase of a string */
char* lower(char* s) {
    
    for (int i = 0; i < strlen(s); i++) {
        s[i] = tolower(s[i]);
    }

    return s;
}

int intArray(int *a, int x, int y) {

    int index = 0;
    for (int i = x; i < y; i++) {
        a[index] = i;
        index++;
    }
    return index;
}

/* [start, end) */
void susbtring(char *res, char* string, int start, int end) {

    char *s = string + start;
    int n = end - start;
    res = strncpy(res, s, n);
    res[n] = '\0';
}

void change(int *a) {
    *a = 10;
}