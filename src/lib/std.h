/* all std-libs, global structs, & variables */

#ifndef STD_H
#define STD_H


#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <ctype.h>
#include <stdbool.h>
#include <math.h>
#include <zlib.h>
#include "kseq.h" /* FASTA-parser */

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