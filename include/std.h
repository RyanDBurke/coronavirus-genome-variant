/*********************************************/
/*                  STD                      */
/* all std-libs, global structs, & variables */
/*********************************************/

#ifndef STD_H
#define STD_H

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <strings.h> /* strcasecmp */
#include <ctype.h>
#include <stdbool.h>
#include <math.h>
#include <zlib.h>
#include "vendor/kseq.h" /* FASTA-parser */

#define FASTA_MAXLINE 512	/* Requires FASTA file lines to be <512 characters */
#define MAXROW 200          /* row bounds on OPT-matrix */
#define MAXCOL 200          /* column bounds on OPT-matrix */

#define FM_SENTINEL '$'    /* end-of-text marker; sorts before every other character */
#define OCC_SAMPLE  64     /* rows between occurrence-count checkpoints */

/* struct for FM-Index over the text T = seq + '$'
   rows 0..length of the (conceptual) burrows-wheeler matrix are the sorted rotations of T;
   row 0 is always the rotation starting at '$' */
typedef struct fm {
    char    *name;          // sequence name
    char    *seq;           // actual sequence string
    int     length;         // length of sequence (T has length + 1 characters)
    int     *suffixArray;   // suffixArray[row] = offset in T where that row's rotation starts
    char    *bwt;           // burrows-wheeler transform: (L)ast column, '$' included
    char    *F;             // (F)irst column of the matrix

    /* fast-rank structures for backward search */
    int     alphabetSize;   // distinct characters in T, '$' included
    int     code[256];      // character -> rank among T's characters (in sort order), -1 if absent
    unsigned char *bwtCodes;// bwt as codes
    int     *C;             // C[c] = # characters in T that sort before c (start of c's block in F)
    int     *occ;           // occ[k * alphabetSize + c] = # of c in bwt[0, k * OCC_SAMPLE)
} FM;

/* struct for an alignment and its score */
typedef struct singleAlignment {
    char    *cigar; // owned; freed by clearAlignments()
    int     score;
    int     pos;   // 1-based reference position (SAM POS)
} Alignment;

/* growable list of alignments */
typedef struct alignments {
    Alignment   *items;
    int         count;
    int         capacity;
} Alignments;

/* tuple for interval [start, end) */
typedef struct interval {
    int start;
    int end; // exclusive
} Interval;

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

#endif
