/**********************************************/
/*                  FM-INDEX                  */
/* associated functions to build the fm-index */
/**********************************************/

#ifndef FMINDEX_H
#define FMINDEX_H

#include "std.h"

/* takes a ref-sequence, builds fmIndex, and writes to output
    * @param fm: fm-index
    * @param reference: .fa file
    * @param output: file output we will write to
 */
int fmIndex(FM *fm, char *reference, char *output);

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

#endif