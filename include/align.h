/***************************************/
/*               ALIGN                 */
/* associated functions to align reads */
/***************************************/

#ifndef ALIGN_H
#define ALIGN_H

#include "std.h"

/* reference bases added on each side of a candidate position before aligning */
#define SLICE_PAD 5

/* read counts from one align() run */
typedef struct alignStats {
    int reads;      // reads in the input
    int mapped;     // reads with at least one alignment
} AlignStats;

/* takes the fmIndex of a ref-sequence, a .fa of single-end reads, and aligns them
    * @param fm: fm-index
    * @param reads: FASTA file of reads (may be gzipped)
    * @param output: SAM file to write, "-" for stdout
    * @param stats: filled with read counts
    * returns EXIT_OK, or the EXIT_* code of the failure (already reported)
 */
int align(FM *fm, const char *reads, const char *output, AlignStats *stats);

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

/* fitting alignment of read at each candidate position; keeps only the best-scoring hits
    * @param hits: best alignments so far for this read (across all seeds); updated in place
    * @param read, readLength: our current read-sequence and its length
    * @param ref, refLength: our reference genome sequence and its length
    * @param refPos: int-array of candidate read start positions (0-based) in the reference
    * @param refPosLength: length of refPos int-array
    * @param gap: our gap penalty
    * @param bestScore: best score so far for this read; updated in place
 */
void alignment(Alignments *hits, char *read, int readLength, char *ref, int refLength, int *refPos, int refPosLength, int gap, int *bestScore);

/* empties an alignment list, freeing its CIGAR strings */
void clearAlignments(Alignments *hits);

/* builds OPT-matrix for a fitting alignment of y inside x, and returns the best score
    (max over the last column, so reference bases before or after the read are free)
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
    * @param offset: set to the slice position (0-based) where the read starts
 */
char *buildCigar(int OPT[MAXROW][MAXCOL], int n, int m, int gap, char *x, char *y, int *offset);

/* return score between two characters: +1 match, -1 mismatch, 0 if either is 'N'
    * @param a: a character in x-axis string (reference genome)
    * @param b: a character in y-axis string (read)
 */
int score(char a, char b);

/* return substring from [start, end)
    * @param res: output for substring
    * @param string: entire string
    * @param start: start of substring
    * @param end: end of string
 */
void substring(char *result, char* string, int start, int end);

/* return run-length encoded CIGAR-string (i.e MMMDDMM -> 3M2D2M)
    * @param ops: one CIGAR op per character, in read order
 */
char *formatCigar(const char *ops);

/* returns our seed skip interval, and intuitively, our seed length */
int seedSkip(int L);

/* min and max */
int min(int a, int b);
int max(int a, int b);
int maxAlign(int a, int b, int c);

/* write to our .SAM file
    * @param sam: our SAM struct
 */
void writeSAM(SAM *sam, FILE *f);

#endif