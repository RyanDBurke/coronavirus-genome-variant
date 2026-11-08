/**************************/
/*         MISC           */
/* miscellanous functions */
/**************************/

#ifndef MISC_H
#define MISC_H

#include "std.h"

/* progress bar */
void pb(bool makeProgressBar, double percentageDone);

/* print commands */
void commands();

/* handle command-line args; returns 0 on success, non-zero for invalid usage
    * @param argc, argv: command-line arguments from main
    * @param ref_seq, indexOut, reads, alignOut: output file paths
    * @param READ: total # of reads (for progress bar), 0 if unknown
    * @param makeProgressBar: whether to draw the progress bar
 */
int handleArgs(int argc, char **argv, char **ref_seq, char **indexOut, char **reads, char **alignOut, double *READ, bool *makeProgressBar);

/* allocation wrappers: print an error and exit if memory runs out, so callers never see NULL */
void *xmalloc(size_t size);
void *xcalloc(size_t count, size_t size);
void *xrealloc(void *ptr, size_t size);

#endif
