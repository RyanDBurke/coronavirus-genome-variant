#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>


#include "fmmap.h"
void test(int x);
void change(int *a);
void susbtring(char *res, char* string, int start, int end);
int intArray(int *a, int x, int y);

int main() {

    for (int i = 0; i < 5; i++) {
        for (int j = 5; j < 10; j++) {
            if ( j == 7) {
                break;
            } else {
                printf("j: %d\n", j);
            }
        }

        printf("i: %d\n", i);
    }

    return 0;

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

int min(int a, int b) {
    return (a > b) ? b : a;
}
int max(int a, int b) {
    return (a > b) ? a : b;
}