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

    /* string path to each relevant file */
    char *ref_seq = argv[1];
    char *indexOut = argv[2];
    char *reads = argv[3];
    char *alignOut = argv[4];

    /* our FM-Index struct */
    FM *fm = malloc(sizeof(FM));

    /* run fmIndex */
    fmIndex(fm, ref_seq, indexOut);

    /* run aligner */
    align(fm, reads, alignOut);

    /* free our fm-index */
    destroy(fm);

   return 0;
}

int fmIndex(FM *fm, char *reference, char *output) {

    FASTAFILE *ffp;
    char *seq;
    char *name;
    int length;

    /* parse .fa file contains our reference sequence and build fm-index */
    ffp = OpenFASTA(reference);
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
        fm->suffixArray = (buildSuffixArray(seq, length));

        /* burrows-wheeler matrix (bwm) */
        fm->bwm = (buildBWM(seq, length));

        /* burrows-wheeler transform (bwt) */
        fm->bwt = malloc(length + 1);
        buildBWT(fm->bwt, fm->bwm, length);

        /* (F)irst and (L)ast columns of bwm */
        fm->F = malloc(length + 1);
        fm->L = malloc(length + 1);
        getFL(fm->F, fm->L, fm->bwm, length);

        /* occTable */
        /* man do I even need this? */
        
        free(seq);
        free(name);
    }

    CloseFASTA(ffp);

    /* write FM-Index to output */
    FILE *f;
    if (length <= 50) {
        f = fopen(output, "w");
        if (f == NULL) {
            printf("error opening file.\n");
            exit(1);
        }
        write(fm, f, length);
    } else {
        printf("> Sequences over 50 in length are not written to file.\n\n");
        printf("> Your sequence length: %d\n\n", length);
        printf("> Use ");
        printf("\033[1;31m");
        printf("./smallFasta.fa ");
        printf("\033[0m");
        printf("to see what a serialized FM-Index looks like.\n");

    }

    /* successful */
    return 0;
}

/* remember, reads have the letter "N" */
int align(FM *fm, char *reads, char *output) {

    FASTAFILE *ffp;
    char *seq;
    char *name;
    int length;

    double ninf = -INFINITY;
    int gap = 5;

    /* parse .fa file containing our n-amount of 100bp reads */
    ffp = OpenFASTA(reads);
    while (ReadFASTA(ffp, &seq, &name, &length)) {

        double best_score = ninf;
        int seedPos = 0;
        int skip = seedSkip(length);
        Alignment A[length]; // is length the correct array size here?

        /* for each (read-length / 5.0) seed (20bp seed in our case, since we have 100bp reads) */
        for (int seedStart = 0; seedStart < length; seedStart += skip) {
            int seedEnd = min(length, seedStart + skip);

            /* finding match interval */
            Interval *interval = malloc(sizeof(Interval)); // I should have to free this, right?
            int matchLength = 0;
            char *seed = malloc(strlen(seq) + 1); // I should have to free this, right?
            substring(seed, seq, seedStart, seedEnd);
            getInterval(interval, &matchLength, fm, seed);

            printf("> %s\n", name);
            printf("seed: %s\n", seed);
            printf("Interval (%d, %d]\n\n", interval->start, interval->end);


            /* free memory */
            free(seed);
            free(interval);
        }

        free(seq);
        free(name);
    }

    CloseFASTA(ffp);

    /* successful */
    return 0;
}

/**************************/
/* AUX FUNTIONS FOR ALIGN */
/**************************/

/* preforms backwards search and return a interval and matchLength */
void getInterval(Interval *interval, int *matchLength, FM *fm, char* seed) {

    /* reset interval to -1 */
    interval->start = -1;
    interval->end = -1;

    /* length of our original reference sequence (coronavirus genome) */
    int refSeqLength = fm->length;

    /* iterate thru suffix array to find matches to our seed */
    for (int i = 0; i < refSeqLength; i++) {

        /* index in coronavirus genome where this particular suffix begins */
        int offset = fm->suffixArray[i];

        /* current coronavirus suffix */
        char *suffix = (fm->seq + offset);

        /* substring of suffix */
        char *subSuffix = malloc(strlen(suffix) + 1); // I should have to free this, right?
        substring(subSuffix, suffix, 0, strlen(seed));

        /* match seed to suffix */
        int match = strcmp(subSuffix, seed);

        /* if they match and this is the first match */
        if (match == 0 && (interval->start == -1)) {
            interval->start = i;
        }

        /* if they don't match, but we've seen a match before -- then this is the end of interval */
        if (match != 0 && (interval->start != -1)) {
            interval->end = i;
            break; // break loop, we have our interval!
        }

        /* free substring created from substring() method */
        free(subSuffix);
    }

    /* matchLength = 0 if no match found, but = strlen(seed) */
    if ((interval->start == -1) && (interval->end == -1)) {
        *matchLength = 0;
    } else {
        *matchLength = strlen(seed);
    }
}

/* return substring of string from [start, end) */
void substring(char* res, char* string, int start, int end) {
    char *temp = string + start;
    int n = end - start;
    res = strncpy(res, temp, n);
    res[n] = '\0';
}

/* seed skip */
int seedSkip(int L) {
    int s = L / 5.0;
    return s;
}

/* min */
int min(int a, int b) {
    return (a > b) ? b : a;
}

/* max */
int max(int a, int b) {
    return (a > b) ? a : b;
}

/*****************************/
/* AUX FUNTIONS FOR FM-INDEX */
/*****************************/

/* builds suffix array */
int *buildSuffixArray(char *seq, int length) {

    /* store suffixes and their offset */
    SA suffixes[length];

    /* starting adding SA objects to array */
    for (int i = 0; i < length; i++) {

        /* offset */
        suffixes[i].offset = i;

        /* build suffix */
        char* currentSuffix = (seq + i);
        suffixes[i].suffix = malloc(strlen(currentSuffix) + 1);
        strcpy(suffixes[i].suffix, currentSuffix);
    }

    /* sort suffixes */
    qsort(suffixes, length, sizeof(SA), &cmpSA);

    /* build suffix array */
    int *sa = malloc(sizeof(int) * length);
    for (int i = 0; i < length; i++) {

        /* store offset in our suffix array */
        sa[i] = suffixes[i].offset;
    }

    return sa;
}

/* build burrows-wheeler matrix */
char **buildBWM(char *seq, int length) {

    /* store rotations and their offset */
    R bwMatrix[length];

    /* starting adding rotations */
    for (int i = 0; i < length; i++) {

        /* offset */
        bwMatrix[i].offset = i;

        /* build rotation */
        char *prefix = malloc(i + 1);
        char *suffix = (seq + i);
        bwMatrix[i].rotation = malloc(length + 1);
        strcpy(bwMatrix[i].rotation, suffix);
        strncpy(prefix, seq, i);
        prefix[i] = '\0'; // add null terminator to prefix (crucial!)
        strcat(bwMatrix[i].rotation, prefix);

        /* free prefix, we don't need it anymore */
        free(prefix);
    }

    /* sort rotations */
    qsort(bwMatrix, length, sizeof(R), &cmpBMW);

    /* stores the actual matrix (rotations sorted) */
    char **matrix = malloc(sizeof(char *) * (length + 1));

    /* fill matrix */
    for (int i = 0; i < length; i++) {
        matrix[i] = malloc(sizeof(char) * length);
        strcpy(matrix[i], bwMatrix[i].rotation);
    }
    
    return matrix;
}

/* builds burrows-wheeler transform (i.e last column) */
void buildBWT(char *bwt, char **BWM, int length) {
    for (int i = 0; i < length; i++) {
        bwt[i] = BWM[i][length - 1];
    }
}

/* gets (F)irst and (L)ast columns of our burrows-wheeler matrix */
void getFL(char *F, char *L, char **BWM, int length) {
    for (int i = 0; i < length; i++) {
        F[i] = BWM[i][0];
        L[i] = BWM[i][length - 1];
    }
}

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

/* write our fm-index to output */
void write(FM *fm, FILE *f, int length) {

    fprintf(f, "************\n");
    fprintf(f, "* FM-INDEX *\n");
    fprintf(f, "************\n\n");

    /* name */
    fprintf(f, "> %s\n", fm->name);

    /* sequence length */
    fprintf(f, "length: %d\n\n", fm->length);

    /* write sequence if its under 512 chars */
    fprintf(f, "Sequence\n");
    fprintf(f, "%s\n\n", fm->seq);

    /* suffix array */
    int lineCount = 0;
    fprintf(f, "Suffix Array\n");
    for (int i = 0; i < length; i++) {
        lineCount++;
        if (i == length - 1) {
            fprintf(f, "%d\n\n", fm->suffixArray[i]);
        } else if (lineCount == 14) {
            fprintf(f, "%d,\n", fm->suffixArray[i]);
            lineCount = 0;
        } else {
            fprintf(f, "%d, ", fm->suffixArray[i]);
        }
    }

    /* bwm */
    fprintf(f, "Burrows-Wheeler Matrix\n");
    for (int i = 0; i < length; i++) {
        if (i == length - 1) {
            fprintf(f, "%s\n\n", fm->bwm[i]);
        } else {
            fprintf(f, "%s\n", fm->bwm[i]);
        }
    }

    /* bwt */
    fprintf(f, "BWT: ");
    fprintf(f, "%s\n\n", fm->bwt);

    /* F and L */
    fprintf(f, "F: ");
    fprintf(f, "%s\n", fm->F);
    fprintf(f, "L: ");
    fprintf(f, "%s\n", fm->L);

    fclose(f);
}

/* prints suffix array */
void printSA(int *sa, int length) {

    printf("Suffix Array\n");
    for (int i = 0; i < length; i++) {
        if (i == length - 1) {
            printf("%d\n\n", sa[i]);
        } else {
            printf("%d, ", sa[i]);
        }
    }
}

/* prints bwm */
void printBWM(char **b, int length) {

    printf("Burrows-Wheeler Matrix\n");
    for (int i = 0; i < length; i++) {
        if (i == length - 1) {
            printf("%s\n\n", b[i]);
        } else {
            printf("%s\n", b[i]);
        }
    }
}

/* prints bwt */
void printBWT(char *bwt, int length) {
    printf("BWT: ");
    printf("%s\n\n", bwt);
}

/* prints F and L */
void printFL(char *F, char *L, int length) {
    printf("F: ");
    printf("%s\n", F);

    printf("L: ");
    printf("%s\n\n", L);
}

/* destroys our FM-Index*/
void destroy(FM *fm) {
    free(fm->bwt);
    free(fm->F);
    free(fm->L);
    free(fm->bwm);
    free(fm->suffixArray);
    free(fm->seq);
    free(fm->name);
    free(fm);
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