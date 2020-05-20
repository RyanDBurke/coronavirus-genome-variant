#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>


#include "fmmap.h"
void test(int x);
void change(int *a);
char *susbtring(char* string, int start, int end);

int main() {

    char *sub = susbtring("ryan", 1, 2);
    printf("%s\n", sub);

    return 0;

}

/* [start, end) */
char *susbtring(char* string, int start, int end) {

    char *s = string + start;
    int n = end - start;
    char *res = malloc(n + 1);
    res = strncpy(res, s, n);
    res[n] = '\0';

    return res;
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