#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>


#include "fmmap.h"
void test(int x);
void change(int *a);
void susbtring(char *res, char* string, int start, int end);
int intArray(int *a, int x, int y);

int main(int argc, char **argv) {

    printf("%s\n", argv[1]);

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