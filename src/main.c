/**************/
/*    MAIN    */
/**************/

#include "align.h"
#include "fmIndex.h"
#include "misc.h"
#include <sys/time.h>
#ifdef _WIN32
#include <fcntl.h>
#include <io.h>
#endif

static double now(void) {
    struct timeval tv;
    gettimeofday(&tv, NULL);
    return tv.tv_sec + tv.tv_usec / 1e6;
}

/* builds the FM-Index of the reference, then aligns the reads to it */
int main(int argc, char **argv) {

#ifdef _WIN32
    /* SAM on stdout must keep "\n" line endings, matching the file output */
    _setmode(_fileno(stdout), _O_BINARY);
#endif

    Options opt;
    switch (parseArgs(argc, argv, &opt)) {
        case ARGS_HELP:
            return EXIT_OK;
        case ARGS_ERROR:
            if (!jsonMode) {
                fprintf(stderr, "run: try 'run --help'\n");
            }
            return EXIT_USAGE;
    }
    quietMode = opt.quiet;
    double start = now();

    /* our FM-Index struct */
    FM *fm = xcalloc(1, sizeof(FM));

    /* run fmIndex, then the aligner */
    AlignStats stats = {0};
    int status = fmIndex(fm, opt.ref, opt.indexOut);
    if (status == EXIT_OK) {
        event("index", "\"length\":%d", fm->length);
        status = align(fm, opt.reads, opt.out, &stats);
    }
    if (status == EXIT_OK) {
        event("done", "\"reads\":%d,\"mapped\":%d,\"unmapped\":%d,\"seconds\":%.2f",
              stats.reads, stats.mapped, stats.reads - stats.mapped, now() - start);
    }

    /* free our fm-index */
    destroy(fm);

    return status;
}
