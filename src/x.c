#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>


#include "fmmap.h"
void test(int x);

int main() {

   double ninf = -INFINITY;
   printf("%f\n", ninf);

    return 0;

}

void test(int x) {
    int h = x / 5.0;
    printf("%d\n", h);
}