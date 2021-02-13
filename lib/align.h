/***************************************/
/*               ALIGN                 */
/* associated functions to align reads */
/***************************************/

#ifndef ALIGN_H
#define ALIGN_H

#include "std.h"

/* takes the fmIndex of a ref-sequence, a .fa of single-end reads, and aligns them 
    * @param fm: fm-index
    * @param reads: FASTA file of 100bp reads
    * @param output: file output we will write to
    * @param READ: total # of reads (for progress bar)
    * @param makePB: to PB or not to PB, that is the question! (PB = progress bar)
 */
int align(FM *fm, char *reads, char *output, double READ, bool makePB);

/**************************/
/* AUX FUNTIONS FOR ALIGN */
/**************************/

/* performs backwards search and returns a interval and matchLength
    * @param interval: interval where seed matched reference genome (query via suffix array)
    * @param matchLength: length of string matched
    * @param fm: fm-index
    * @param partialSeq: substring of our current read sequence from seedStart:seedEnd
 */
void getInterval(Interval *interval, int *matchLength, FM *fm, char* seed);

/* returns refPos-array length and reference positions to its parameter
    * @param refPos: int-array of all reference positions
    * @param interval: interval where seed matched reference genome
    * @param matchLength: length of string matched to seed
    * @param fm: fm-index
    * @param seedEnd: index where current seed ends
 */
int referencePos(int *refPos, Interval *interval, int matchLength, FM *fm, int seedEnd);

/* return a single alignment struct to output parameter A, and length of alignments-array
    * @param alignments: an array of Alignment structs
    * @param read: our current read-sequence
    * @param ref: our reference genome sequence
    * @param refPos: int-array of all reference positions in our reference genome
    * @param refPosLength: length of refPos int-array
    * @param gap: our gap penalty
    * @param bestScore: our current best score
 */ 
int alignment(Alignment alignments[], char *read, char *ref, int *refPos, int refPosLength, int gap, double bestScore);

/* builds OPT-matrix, and returns OPT[n][m] -- which is the edit distance between x and y
    * @param OPT: our OPT-matrix
    * @param x: our slice from the reference genome
    * @param y: our read
    * @param n, m: lengths of strings x and y. Intuitively, also our x/y-axis
    * @param gap: our gap penalty
 */
int editDistance(int OPT[MAXROW][MAXCOL], char *x, char *y, int n, int m, int gap);

/* returns a CIGAR-string of our alignment
    
    * @param OPT: our OPT-matrix
    * @param n, m: lengths of strings x and y. Intuitively, also our x/y-axis
    * @param gap: our gap penalty
    * @param x: our slice from the reference genome
    * @param y: our read
 */
char *buildCigar(int OPT[MAXROW][MAXCOL], int n, int m, int gap, char *x, char *y, int *offset);

/* return score between two characters
    * @param a: a character in x-axis string (reference genome)
    * @param b: a character in y-axis string (read)
    * @param gap: our gap penalty
 */
int score(char a, char b, int gap);

/* return substring from [start, end)
    * @param res: output for substring
    * @param string: entire string
    * @param start: start of substring
    * @param end: end of string
 */
void substring(char *result, char* string, int start, int end);

/* return correctly formatted CIGAR-string (i.e 3M2D2M)
    * @param cigar: incorrectly formatted cigar-string (i.e MMMDDMM)
    * @param length: length of cigar-string passed in
 */
char *formatCigar(char *cigar, int length);

/* returns our seed skip interval, and intuitively, our seed length */
int seedSkip(int L);

/* min and max */
int min(int a, int b);
int max(int a, int b);
int maxAlign(int a, int b, int c);
int minAlign(int a, int b, int c);

/* print methods 
    * @param OPT: our OPT-matrix
    * @param n, m: our genome-string lengths + 1
    * @param x, y: our reference genome and read genome, respectively
 */
void printMatrix(int OPT[MAXROW][MAXCOL], int n, int m, char *x, char *y);

/* write to our .SAM file
    * @param sam: our SAM struct
 */
void writeSAM(SAM *sam, FILE *f);

/* free our SAM struct */
void destroySAM(SAM *sam);

/* destroy our alignment array */
void destroyAlignment(Alignment *A);

#endif