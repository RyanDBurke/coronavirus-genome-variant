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
int max(int a, int b);
char *formatCigar(char *cigar, int length);

int main(int argc, char **argv) {

    char *x = "AAGGTATGAATC"; // 12
    char *y = "AACGTTGAC";   // 9
    int n = strlen(x);      // 12
    int m = strlen(y);      // 9
    int matrix[MAXROW][MAXCOL];
    int gap = -3;
    
    int edit = editDistance(matrix, x, y, n, m, gap);

    // printMatrix(matrix, n, m, x, y);

    char *z = "MMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMM";

    int offset = -1;
    char *cigar = buildCigar(matrix, n, m, gap, x, y, &offset);
    char *trueCigar = formatCigar(z, strlen(z));
    printf("CIGAR: \t%s\n", z);
    printf("TRUE CIGAR: \t%s\n", trueCigar);
    
    free(trueCigar);
    free(cigar);

    char array[64];
    int myInteger = 10;
    sprintf( array, "%d", myInteger );

    return 0;

}

char *formatCigar(char *cigar, int length) {

    /* just make sure cigar input is uppercase */
    // cigar = upper(cigar);

    char *result = malloc(length + 1);
    int resultIndex = 0;

    int count = 0;
    char currentChar = '%';
    for (int i = 0; i < length; i++) {
        currentChar = cigar[i];

        for (int j = i; j < length; j++) {
            if (cigar[j] != currentChar) {
                if (count > 9) {
                    char array[64];
                    int myInteger = count;
                    sprintf(array, "%d", myInteger );

                    for (int k = 0; k < strlen(array); k++) {
                        result[resultIndex] = array[k];
                        resultIndex++;
                    }

                    // strcat(result, array);
                    // resultIndex += strlen(array);
                } else {
                    result[resultIndex] = count + '0';
                    resultIndex++;
                }
                
                result[resultIndex] = currentChar;
                resultIndex++;

                if (cigar[j + 1] == '\0') {
                    j = length;
                }
                break;
            } else {
                count++;

                if (cigar[j + 1] == '\0') {
                    if (count > 9) {
                        char array[64];
                        int myInteger = count;
                        sprintf( array, "%d", myInteger );

                        for (int k = 0; k < strlen(array); k++) {
                            result[resultIndex] = array[k];
                            resultIndex++;
                        }

                        // strcat(result, array);
                        // resultIndex += strlen(array);
                    } else {
                        result[resultIndex] = count + '0';
                        resultIndex++;
                    }

                    result[resultIndex] = currentChar;
                    resultIndex++;
                    j = length;
                }
            }
            i = j;
        }



        count = 0;
    }

    result[resultIndex] = '\0';
    printf("cigar: %s\n", result);

    return result;
}

int max(int a, int b) {
    return (a > b) ? a : b;
}

/* returns CIGAR string  */
char *buildCigar(int OPT[MAXROW][MAXCOL], int n, int m, int gap, char *x, char *y, int *offset) {

    /* find out where we start our backtrace, and set n */
    int temp = n;
    int maxScore = -1;
    for (int i = 0; i <= temp; i++) {
        if (OPT[i][m] >= maxScore) {
            maxScore = OPT[i][m];
            n = i;
        }
    }

    /* result array we'll return */
    int traceIndex = 0;


    char *cigar = malloc(n * m);
    while (true) {

        /* we've reached matrix[0][0] */
        /* for fitting alignment it ends when we reach row = 0 */
        if (m == 0) {
            *offset = n;
            break;
        }

        int up = gap + OPT[n - 1][m];
        int diagonal = score(x[n - 1], y[m - 1], gap) + OPT[n - 1][m - 1];
        int left = gap + OPT[n][m - 1];

        int maxA = maxAlign(up, diagonal, left);
        if (maxA == diagonal) {
            n = n - 1;
            m = m - 1;
            cigar[traceIndex] = 'M';
            traceIndex++;
        } else if (maxA == left) {
            m = m - 1;
            cigar[traceIndex] = 'I';
            traceIndex++;
        } else {
            n = max(n - 1, 0);
            cigar[traceIndex] = 'D';
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
        return 1;
    } else if (a != b) {
        return -1;
    } else if (a == '-'  || b == '-') {
        return gap;
    } else {
        return 0;
    }
}

void printMatrix(int matrix[MAXROW][MAXCOL], int n, int m, char *x, char *y) {

    for (int i = 0; i <= m + 1; i++) {
        if (i == 0) {
            printf("\t");
        } else if (i == 1) {
            printf("Y-axis\t");
        } else {
            printf("%c\t", y[i - 2]);
        }
    }

    printf("\n");

    for (int i = 0; i <= n; i++) {     
        if (i > 0) {
            printf("%c\t", x[i - 1]);
        } else {
            printf("X-axis\t");
        }
        for (int j = 0; j <= m; j++) {
            printf("%d\t", matrix[i][j]);
        }
        printf("\n");
    }
}

/* builds OPT-matrix, and returns OPT[n][m] -- which is the edit distance between x and y */
int editDistance(int OPT[MAXROW][MAXCOL], char *x, char *y, int n, int m, int gap) {

    /* initial [0][0] is always 0 */
    OPT[0][0] = 0;

    /* add our initial gap penalties to the first column of each row */
    /* cost-free ends */
    for (int i = 1; i <= n; i++) {
        OPT[i][0] = 0;
    }

    /* add our initial gap penalties to the first row of each column */
    for (int j = 1; j <= m; j++) {
        OPT[0][j] = j * gap;
    }

    /* fill matrix */
    /* whether i make this <= or < it seems to work? Why? */
    for (int i = 1; i <= n; i++) {
        for (int j = 1; j <= m; j++) {

            char a = x[i - 1];
            char b = y[j - 1];

            /* calculate score */
            OPT[i][j] = maxAlign(
                    (score(a, b, gap) + OPT[i - 1][j - 1]),
                    (gap + OPT[i - 1][j]),
                    (gap + OPT[i][j - 1])
                );

            /*
            printf("Comparing: %c -> %c\n", a, b);
            printf("Score 1: %d\n", (score(a, b, gap) + OPT[i - 1][j - 1]));
            printf("Score 2: %d\n", (gap + OPT[i - 1][j]));
            printf("Score 3: %d\n", (gap + OPT[i][j - 1]));
            printf("Final Score: %d\n\n", OPT[i][j]);
            */
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