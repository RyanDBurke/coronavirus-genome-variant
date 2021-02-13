/**************/
/*    MAIN    */
/**************/

#include "../lib/align.h"
#include "../lib/fmIndex.h"
#include "../lib/misc.h"

/* keep a global note on # of reads */
double READ = 0;
bool makeProgressBar = false;

/* runs align and fmIndex */
int main(int argc, char **argv) {

    /* string path to each relevant file */
    char *ref_seq = argv[1];
    char *indexOut = argv[2];
    char *reads = argv[3];
    char *alignOut = argv[4];

    /* Valid commands */
    commands();

    /* handle args */
    handleArgs(argv, &ref_seq, &indexOut, &reads, &alignOut, &READ, &makeProgressBar);

    /* our FM-Index struct */
    FM *fm = malloc(sizeof(FM));

    /* run fmIndex */
    fmIndex(fm, ref_seq, indexOut);

    /* run aligner */
    align(fm, reads, alignOut, READ, makeProgressBar);

    /* free our fm-index */
    destroy(fm);

   return 0;
}
