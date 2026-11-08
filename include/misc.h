/****************************************************/
/*                     MISC                         */
/* command-line handling, messages, and allocation */
/****************************************************/

#ifndef MISC_H
#define MISC_H

#include "std.h"

/* everything main() needs to know from the command line */
typedef struct options {
    const char  *ref;           // reference FASTA
    const char  *reads;         // reads FASTA
    const char  *out;           // SAM output path, "-" for stdout
    const char  *indexOut;      // FM-Index text dump (references up to 50 bp)
    double      readCount;      // total reads for the progress bar, 0 if unknown
    bool        quiet;          // no progress bar or status messages
} Options;

/* parseArgs results */
#define ARGS_OK    0
#define ARGS_ERROR 1
#define ARGS_HELP  2

/* parse the command line into opt; prints its own error messages
    * returns ARGS_OK, ARGS_ERROR (invalid usage) or ARGS_HELP (--help was printed)
 */
int parseArgs(int argc, char **argv, Options *opt);

/* print usage to f */
void usage(FILE *f);

/* status messages go to stderr (stdout may be carrying SAM), and are silenced by --quiet */
extern bool quietMode;
void info(const char *fmt, ...) __attribute__((format(printf, 1, 2)));

/* highlight escape codes, or "" when stderr is not a terminal */
const char *highlight(void);
const char *highlightEnd(void);

/* progress bar on stderr: redraws only when the whole-number percentage changes */
void pb(bool makeProgressBar, double percentageDone);

/* allocation wrappers: print an error and exit if memory runs out, so callers never see NULL */
void *xmalloc(size_t size);
void *xcalloc(size_t count, size_t size);
void *xrealloc(void *ptr, size_t size);

#endif
