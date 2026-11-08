/**************/
/*    MAIN    */
/**************/

#include "align.h"
#include "fmIndex.h"
#include "misc.h"

/* builds the FM-Index of the reference, then aligns the reads to it */
int main(int argc, char **argv) {

    Options opt;
    switch (parseArgs(argc, argv, &opt)) {
        case ARGS_HELP:
            return 0;
        case ARGS_ERROR:
            fprintf(stderr, "run: try 'run --help'\n");
            return 1;
    }
    quietMode = opt.quiet;

    /* our FM-Index struct */
    FM *fm = xcalloc(1, sizeof(FM));

    /* run fmIndex, then the aligner (progress bar only when the read count is known) */
    int status = fmIndex(fm, opt.ref, opt.indexOut);
    if (status == 0) {
        status = align(fm, opt.reads, opt.out, opt.readCount, opt.readCount > 0);
    }

    /* free our fm-index */
    destroy(fm);

    return status == 0 ? 0 : 1;
}
