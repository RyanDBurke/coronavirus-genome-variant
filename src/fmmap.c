#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <ctype.h>
#include <stdbool.h>

#include "fmmap.h"

/* runs align and fmIndex */
/* ./fmmap ref_seq indexOutput reads alignOutput */
int main(int argc, char **argv) {

    /* invalid number of arguments
    if (argc != 5) {
        printf("invalid number of arguments.");
        return 0;
    }
    */

    /* string path to each file */
    char *ref_seq = argv[1];
    char *indexOut = argv[2];
    // char *reads = argv[3];
    // char *alignOut = argv[4];

    /* fmIndex */
    fmIndex(ref_seq, indexOut);

   return 0;
}

int fmIndex(char *reference, char *output) {

    FASTAFILE *ffp;
    char *seq;
    char *name;
    int length;

    /* parse .fa file contains our reference sequence and build fm-index */
    ffp = OpenFASTA(reference);
    FM *fm = malloc(sizeof(FM));
    while (ReadFASTA(ffp, &seq, &name, &length)) {

        /* name */
        fm->name = malloc(strlen(name) + 1);
        strcpy(fm->name, name);

        /* seq */
        fm->seq = malloc(strlen(seq) + 1);
        strcpy(fm->seq, seq);
        
        /* length of seq */
        fm->length = length;

        /* suffix array */
        fm->suffixArray = (buildSuffixArray(seq, length + 1));
        // printSA(fm->suffixArray, length + 1);

        /* burrows-wheeler matrix (bwm) */
        fm->bwm = (bw(seq, length));
        // printBWM(fm->bwm, length);

        /* occTable */

        free(seq);
        free(name);
    }

    CloseFASTA(ffp);

    /* free our fm-index */
    free(fm->bwm);
    free(fm->suffixArray);
    free(fm->seq);
    free(fm->name);
    free(fm);

    return 0;
}

int align() {
    return 0;
}

/****************/
/* AUX FUNTIONS */
/****************/


/* comparison sort function for suffix arrays */
int cmpSA(const void *a, const void *b) {
    const SA *da = (const SA *) a;
    const SA *db = (const SA *) b;

    return strcmp(da->suffix, db->suffix);
}

/* comparison sort function for bwm */
int cmpBMW(const void *a, const void *b) {
    const R *da = (const R *) a;
    const R *db = (const R *) b;

    return strcmp(da->rotation, db->rotation);
}

/* builds suffix array */
int *buildSuffixArray(char *seq, int length) {

    /* store suffixes and their offset */
    SA suffixes[length];

    /* starting adding SA objects to array */
    for (int i = 0; i < length; i++) {

        /* offset */
        suffixes[i].offset = i;

        /* build suffix */
        suffixes[i].suffix = malloc(strlen(seq + i) + 1);
        strcpy(suffixes[i].suffix, (seq + i));

        /* add null character */
        /* do i need this */
        (suffixes[i].suffix)[length] = '\0';
    }

    /* sort suffixes */
    qsort(suffixes, length, sizeof(SA), &cmpSA);

    /* build suffix array */
    int *sa = malloc(sizeof(int) * length);
    for (int i = 0; i < length; i++) {
        sa[i] = suffixes[i].offset;
        free(suffixes[i].suffix);
    }

    return sa;
}

/* build burrows-wheeler matrix */
char **bw(char *seq, int length) {

    /* store rotations and their offset */
    R bwMatrix[length];

    /* starting adding rotations */
    for (int i = 0; i < length; i++) {

        /* offset */
        bwMatrix[i].offset = i;

        /* build rotation */
        char *prefix = malloc(i + 1);
        bwMatrix[i].rotation = malloc(length);
        strcpy(bwMatrix[i].rotation, (seq + i));
        strncpy(prefix, seq, i);
        strcat(bwMatrix[i].rotation, prefix);

        /* remember to add null character! */
        (bwMatrix[i].rotation)[length] = '\0';

        /* free used memory */
        free(prefix);
    }

    /* sort rotations */
    qsort(bwMatrix, length, sizeof(R), &cmpBMW);

    /* stores the actual matrix, with just the rotations */
    char **matrix = malloc(sizeof(char *) * length);

    /* fill matrix */
    for (int i = 0; i < length; i++) {
        matrix[i] = malloc(sizeof(char) * length);
        // *matrix = malloc(sizeof(char) * length);
        strcpy(matrix[i], bwMatrix[i].rotation);
        free(bwMatrix[i].rotation);
    }
    
    return matrix;
}

/* prints suffix array */
void printSA(int *sa, int length) {
    for (int i = 0; i < length; i++) {
        printf("%d\n", sa[i]);
    }
}

/* prints bwm */
void printBWM(char **b, int length) {
    for (int i = 0; i < length; i++) {
        printf("%s\n", b[i]);
    }
}


/****************/
/* FASTA PARSER */
/****************/

FASTAFILE *
OpenFASTA(char *seqfile)
{
  FASTAFILE *ffp;

  ffp = malloc(sizeof(FASTAFILE));
  ffp->fp = fopen(seqfile, "r");              /* Assume seqfile exists & readable!   */
  if (ffp->fp == NULL) { free(ffp); return NULL; } 
  if ((fgets(ffp->buffer, FASTA_MAXLINE, ffp->fp)) == NULL)
    { free(ffp); return NULL; }
  return ffp;
}

int
ReadFASTA(FASTAFILE *ffp, char **ret_seq, char **ret_name, int *ret_L)
{
  char *s;
  char *name;
  char *seq;
  int   n;
  int   nalloc;
  
  /* Peek at the lookahead buffer; see if it appears to be a valid FASTA descline.
   */
  if (ffp->buffer[0] != '>') return 0;    

  /* Parse out the name: the first non-newline token after the >
   */
  s  = strtok(ffp->buffer+1, "\n");
  name = malloc(sizeof(char) * (strlen(s)+1));
  strcpy(name, s);

  /* Everything else 'til the next descline is the sequence.
   * Note the idiom for dynamic reallocation of seq as we
   * read more characters, so we don't have to assume a maximum
   * sequence length.
   */
  seq = malloc(sizeof(char) * 128);     /* allocate seq in blocks of 128 residues */
  nalloc = 128;
  n = 0;
  while (fgets(ffp->buffer, FASTA_MAXLINE, ffp->fp))
    {
      if (ffp->buffer[0] == '>') break;	/* a-ha, we've reached the next descline */

      for (s = ffp->buffer; *s != '\0'; s++)
	{
	  if (! isalpha(*s)) continue;  /* accept any alphabetic character */

	  seq[n] = *s;                  /* store the character, bump length n */
	  n++;
	  if (nalloc == n)	        /* are we out of room in seq? if so, expand */
	    {			        /* (remember, need space for the final '\0')*/
	      nalloc += 128;
	      seq = realloc(seq, sizeof(char) * nalloc);
	    }
	}
    }
  seq[n] = '\0';

  *ret_name = name;
  *ret_seq  = seq;
  *ret_L    = n;
  return 1;
}      

void
CloseFASTA(FASTAFILE *ffp)
{
  fclose(ffp->fp);
  free(ffp);
}