/**************/
/*    MAIN    */
/**************/

#include "../lib/align.h"
#include "../lib/fmIndex.h"
#include "../lib/misc.h"

/* runs align and fmIndex */
int main(int argc, char **argv) {

    /* string path to each relevant file */
    char *ref_seq, *indexOut, *reads, *alignOut;
    double READ = 0;
    bool makeProgressBar = false;

    /* handle args */
    if (handleArgs(argc, argv, &ref_seq, &indexOut, &reads, &alignOut, &READ, &makeProgressBar) != 0) {
        fprintf(stderr, "invalid arguments\n");
        commands();
        return 1;
    }

    /* our FM-Index struct */
    FM *fm = calloc(1, sizeof(FM));

    /* run fmIndex, then the aligner */
    int status = fmIndex(fm, ref_seq, indexOut);
    if (status == 0) {
        status = align(fm, reads, alignOut, READ, makeProgressBar);
    }

    /* free our fm-index */
    destroy(fm);

    return status == 0 ? 0 : 1;
}
