/* fmmap.h */

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <ctype.h>
#include <stdbool.h>

#define FASTA_MAXLINE 512 /* Requires FASTA file lines to be <512 characters */

/* struct for FM-Index */
typedef struct fm {
    char *name;
    char *seq;
    int length;
    int *suffixArray;
    char **bwm;
    char *F;
    char *L;
    // occTable
} FM;

/* suffix array struct */
typedef struct suffixArray {
    int offset;
    char *suffix;
} SA;

/* rotation struct for burrows-wheeler matrix */
typedef struct Rotation {
    int offset;
    char *rotation;
} R;

/* FASTA files */
typedef struct fastafile_s {
  FILE *fp;
  char  buffer[FASTA_MAXLINE];
} FASTAFILE;

/* FASTA parser */
extern FASTAFILE *OpenFASTA(char *seqfile);
extern int        ReadFASTA(FASTAFILE *fp, char **ret_seq, char **ret_name, int *ret_L);
extern void       CloseFASTA(FASTAFILE *ffp);

/* takes a ref-sequence, builds fmIndex, and writes to output
    * @param reference: .fa file
    * @param output: file output we will write to
 */
int fmIndex(char *reference, char *output);

/* takes the fmIndex of a ref-sequence, a .fa of reads, and aligns them */
int align();

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

/* builds occurence table
    * @param: occF: address to a map-like struct <character: occurence #>
    * @param: occL: address to a map-like struct <character: occurence #>
    * @param: BWM: our burrows-wheeler matrix (useful for extracting (F)irst and (L)ast columns)
    * pass in &occF and &occL
 */
void buildOccTable(char **occF, char **occL, char **BWM); // char for now, change later

/* gets (F)irst and (L)ast columns of our burrows-wheeler matrix
    * @param F: char array of the first column
    * @param L: char array of last column
    * @param BWM: our burrows-wheeler matrix
    * we'd pass in &F and &L
 */
void getFL(char **F, char **L, char **BWM, int length);

/* comparison sort functions for suffix arrays and BWM */
int cmpSA(const void *a, const void *b);
int cmpBMW(const void *a, const void *b);

/* print methods */
void printSA(int *sa, int length);
void printBWM(char **b, int length);
void printFL(char *F, char *L, int length);