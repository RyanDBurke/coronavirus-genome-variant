/* fmmap.h */

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <ctype.h>
#include <stdbool.h>
#include <math.h>
#include <zlib.h>

#define FASTA_MAXLINE 512	/* Requires FASTA file lines to be <512 characters */
#define MAXROW 200
#define MAXCOL 200

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

    /* not implemented yet */
    int     *F_rank;       // fast-rank calculation for (F)irst column
    int     *L_rank;       // fast-rank calculation for (L)ast column
} FM;

/* struct for an alignment and its score */
typedef struct singleAlignment {
    char    *cigar; // needs to be free'd in each alignment
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

/* SAM format */
typedef struct sam {
    char    *QNAME;
    int     FLAG;
    char    *RNAME;
    int     POS;
    int     MAPQ;
    char    *CIGAR;
    char    *RNEXT;
    int     PNEXT;
    int     TLEN;
    char    *SEQ;
    char    *QUAL;
} SAM;

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
void writeFM(FM *fm, FILE *f, int length);

/* print methods */
void printSA(int *sa, int length);
void printBWM(char **b, int length);
void printBWT(char *bwt, int length);
void printFL(char *F, char *L, int length);

/* destroys our FM-Index */
void destroy(FM *fm);

/********/
/* MISC */
/********/

/* string to lowercase */
char* lower(char* s);

/* string to uppercase*/
char* upper(char* s);

/* progress bar */
void pb(bool makeProgressBar, double percentageDone);