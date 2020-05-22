#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <ctype.h>
#include <stdbool.h>
#include <math.h>

const int row = 200;
const int col = 200;

// #include "fmmap.h"
void test(int x);
void change(int *a);
void susbtring(char *res, char* string, int start, int end);
int intArray(int *a, int x, int y);
int maxAlign(int a, int b, int c);
int score(char a, char b, int gap);
int editDistance(int OPT[row][col], char *x, char *y, int n, int m, int gap);
int minAlign(int a, int b, int c);

int main(int argc, char **argv) {

    char *x = "AAGGTATGAATC";
    char *y = "AACGTTGAC";
    int n = strlen(x) + 1;
    int m = strlen(y) + 1;
    int matrix[row][col];
    int gap = 2;

    /* our edit-distance */
    int edit = editDistance(matrix, x, y, n, m, gap);

    printf("edit distance: %d\n", edit);

    for (int i = 0; i < n; i++) {
        for (int j = 0; j < m; j++) {
            printf("%d\t", matrix[i][j]);
        }
        printf("\n");
    }


    return 0;

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

/* builds OPT-matrix, and returns OPT[n][m] -- which is the edit distance between x and y */
int editDistance(int OPT[row][col], char *x, char *y, int n, int m, int gap) {

    /* add our initial gap penalties to the first column of each row */
    for (int i = 0; i < n; i++) {
        OPT[i][0] = 0;
    }

    /* add our initial gap penalties to the first row of each column */
    for (int j = 0; j < m; j++) {
        OPT[0][j] = 0;
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

void test(int x) {
    int h = x / 5.0;
    printf("%d\n", h);
}