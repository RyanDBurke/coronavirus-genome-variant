#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>


#include "fmmap.h"
void test(int x);
void change(int *a);
void susbtring(char *res, char* string, int start, int end);

int main() {

    char *string = "ryanryanryanryan";
    char *sub = malloc(strlen(string) + 1);
    susbtring(sub, string, 1, 7);
    printf("%s\n", sub);
    int l = strlen(sub);
    printf("%d\n", l);
    

    free(sub);

    return 0;

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