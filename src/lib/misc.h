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

#endif