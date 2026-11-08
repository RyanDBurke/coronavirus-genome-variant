/**********************************************/
/*                  FM-INDEX                  */
/* associated functions to build the fm-index */
/**********************************************/

#ifndef FMINDEX_H
#define FMINDEX_H

#include "std.h"

/* takes a ref-sequence, builds fmIndex, and writes it to output (for short references)
    * @param fm: fm-index
    * @param reference: .fa file
    * @param output: file output we will write to
 */
int fmIndex(FM *fm, char *reference, char *output);

/*****************************/
/* AUX FUNTIONS FOR FM-INDEX */
/*****************************/

/* builds the suffix array, BWT, F column, and rank tables for fm->seq / fm->length
    * @param fm: fm-index with seq and length already set
 */
void buildIndex(FM *fm);

/* builds suffix array of seq + '$'
    * @param seq: our sequence
    * @param length: length of seq
    * returns length + 1 offsets in sorted-suffix order; entry 0 is always length ('$')
 */
int *buildSuffixArray(char *seq, int length);

/* number of occurrences of character code c in bwt[0, row) (fast rank)
    * @param fm: fm-index
    * @param c: character code (fm->code[character])
    * @param row: exclusive end row
 */
int occCount(FM *fm, int c, int row);

/* FM-Index backward search for an exact pattern
    * @param fm: fm-index
    * @param pattern, m: the pattern and its length
    * @param interval: set to the suffix-array rows [start, end) that begin with pattern
    * returns true if the pattern occurs
 */
bool backwardSearch(FM *fm, const char *pattern, int m, Interval *interval);

/* comparison sort function for suffix array offsets */
int cmpSA(const void *a, const void *b);

/* write our fm-index to output */
void writeFM(FM *fm, FILE *f);

/* frees the index structures built by buildIndex (not name/seq) */
void freeIndex(FM *fm);

/* destroys our FM-Index */
void destroy(FM *fm);

#endif
