/**************************/
/*         MISC           */
/* miscellanous functions */
/**************************/

#ifndef MISC_H
#define MISC_H

#include "std.h"

/* string to lowercase */
char* lower(char* s);

/* string to uppercase*/
char* upper(char* s);

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

#endif
