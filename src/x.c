#include <stdlib.h>
#include <stdio.h>
#include <string.h>


#include "fmmap.h"

void test(char **f);

int main() {

    FM *fm = malloc(sizeof(FM));
    fm->F = malloc(6);
    strcpy(fm->F, "hello");
    // printf("%s\n", fm->F);
    test(&(fm->F));
    // printf("%s\n", fm->F);
    free(fm->F);
    free(fm); // why not free?

    return 0;

}

void test(char **F) {
    *F = "ryan";
}