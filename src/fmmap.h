/* fmmap.h */

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <ctype.h>
#include <stdbool.h>
#include <math.h>

#define FASTA_MAXLINE 512	/* Requires FASTA file lines to be <512 characters */

/* struct for FM-Index */
typedef struct fm {
    char    *name;         // sequence name
    char    *seq;          // actual sequence string
    int     length;        // length of sequence
    int     *suffixArray;  // suffix array of sequence
    char    **bwm;         // burrows-wheeler matrix
    char    *bwt;          // burrows-wheeler transform
    char    *F;            // (F)irst column of bwm
    char    *L;            // (L)ast column of bwm
    // occTable
} FM;

/* struct for an alignment and its score */
typedef struct singleAlignment {
    char    *alignmentX; // because we would need two, correct?
    char    *alignmentY;
    int     score;
    int     pos;
} Alignment;

/* tuple for interval [start, end) */
typedef struct interval {
    int start;
    int end; // exclusive
} Interval;

/* suffix array struct */
typedef struct suffixArray {
    int     offset;
    char    *suffix;
} SA;

/* rotation struct for burrows-wheeler matrix */
typedef struct Rotation {
    int     offset;
    char    *rotation;
} R;

/* FASTA files */
typedef struct fastafile_s {
  FILE  *fp;
  char  buffer[FASTA_MAXLINE];
} FASTAFILE;

/******************/
/* MAIN FUNCTIONS */
/******************/

/* takes a ref-sequence, builds fmIndex, and writes to output
    * @param fm: fm-index
    * @param reference: .fa file
    * @param output: file output we will write to
 */
int fmIndex(FM *fm, char *reference, char *output);

/* takes the fmIndex of a ref-sequence, a .fa of single-end reads, and aligns them 
    * @param fm: fm-index
    * @param reads: FASTA file of 100bp reads
    * @param output: file output we will write to
 */
int align(FM *fm, char *reads, char *output);

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

/* return a single alignment struct to output parameter A
    * @param A: a single alignment struct, containing updated values
    * @param read: our current read-sequence
    * @param ref: our reference genome sequence
    * @param refPos: int-array of all reference positions in our reference genome
    * @param gap: our gap penalty
 */ 
void alignment(Alignment *A, char *read, char *ref, int *refPos, int gap);

/* return substring from [start, end)
    * @param res: output for substring
    * @param string: entire string
    * @param start: start of substring
    * @param end: end of string
 */
void substring(char *result, char* string, int start, int end);

/* returns our seed skip interval, and intuitively, our seed length */
int seedSkip(int L);

/* min and max */
int min(int a, int b);
int max(int a, int b);

/*****************************/
/* AUX FUNTIONS FOR FM-INDEX */
/*****************************/

/* builds suffix array
    * @param seq: .fa file
    * @param length: length of our seq (.fa file)
 */
int *buildSuffixArray(char *seq, int length);

/* builds burrows-wheeler matrix
    * @param seq: .fa file
    * @param length: length of our seq (.fa file)
 */
char **buildBWM(char *seq, int length);

/* builds burrows-wheeler transform (i.e last column)
    * @param BWM: our burrows-wheeler matrix
    * @param length: length of sequence
 */
void buildBWT(char *bwt, char **BWM, int length);

/* builds occurence table
    * @param: occF: address to a map-like struct <character: occurence #>
    * @param: occL: address to a map-like struct <character: occurence #>
    * @param: BWM: our burrows-wheeler matrix (useful for extracting (F)irst and (L)ast columns)
 */
void buildOccTable(char **occF, char **occL, char **BWM); // char for now, change later

/* return (F)irst and (L)ast columns of our burrows-wheeler matrix in F and L
    * @param F: char array of the first column
    * @param L: char array of last column
    * @param BWM: our burrows-wheeler matrix
    * @param length: length of sequence
 */
void getFL(char *F, char *L, char **BWM, int length);

/* comparison sort functions for suffix arrays and BWM */
int cmpSA(const void *a, const void *b);
int cmpBMW(const void *a, const void *b);

/* write our fm-index to output */
void write(FM *fm, FILE *f, int length);

/* print methods */
void printSA(int *sa, int length);
void printBWM(char **b, int length);
void printBWT(char *bwt, int length);
void printFL(char *F, char *L, int length);

/* destroys our FM-Index */
void destroy(FM *fm);

/****************/
/* FASTA PARSER */
/****************/

extern FASTAFILE *OpenFASTA(char *seqfile);
extern int        ReadFASTA(FASTAFILE *fp, char **ret_seq, char **ret_name, int *ret_L);
extern void       CloseFASTA(FASTAFILE *ffp);