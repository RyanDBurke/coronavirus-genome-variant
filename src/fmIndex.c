#include "fmIndex.h"
#include "misc.h"
KSEQ_INIT(gzFile, gzread)

/* references up to this length get their full FM-Index written to a text file */
#define FM_WRITE_LIMIT 50

int fmIndex(FM *fm, const char *reference, const char *output) {

    /* fasta stuff */
    if (reference == NULL) {
        fprintf(stderr, "no reference file given\n");
        return 1;
    }
    gzFile fp = gzopen(reference, "r");
    if (fp == NULL) {
        fprintf(stderr, "error opening reference file: %s\n", reference);
        return 1;
    }
    kseq_t *seq = kseq_init(fp);

    /* parse .fa file and build the fm-index from its first sequence */
    if (kseq_read(seq) < 0 || seq->seq.l == 0) {
        fprintf(stderr, "no sequence found in reference file: %s\n", reference);
        kseq_destroy(seq);
        gzclose(fp);
        return 1;
    }

    /* name (SAM RNAME), taken from the FASTA header; fall back if the header has none */
    const char *name = (seq->name.l > 0) ? seq->name.s : "ref";
    fm->name = xmalloc(strlen(name) + 1);
    strcpy(fm->name, name);

    /* seq */
    fm->length = (int) seq->seq.l;
    fm->seq = xmalloc(fm->length + 1);
    memcpy(fm->seq, seq->seq.s, fm->length + 1);

    kseq_destroy(seq);
    gzclose(fp);

    /* suffix array, burrows-wheeler transform, and rank tables */
    buildIndex(fm);

    /* write FM-Index to output */
    if (fm->length <= FM_WRITE_LIMIT) {
        FILE *f = fopen(output, "w");
        if (f == NULL) {
            fprintf(stderr, "error opening FM-Index output file: %s\n", output);
            return 1;
        }
        writeFM(fm, f);
        fclose(f);
        info("> FM-Index for \"%s\" written to %s%s%s\n", reference, highlight(), output, highlightEnd());
    } else {
        info("> Reference is %d bp; the FM-Index is only written to a file for references up to %d bp (try ./run)\n",
             fm->length, FM_WRITE_LIMIT);
    }

    return 0;
}

/*****************************/
/* AUX FUNTIONS FOR FM-INDEX */
/*****************************/

/* text being sorted by cmpSA (qsort has no context argument) */
static const char *sortText;

/* comparison sort function for suffix array offsets: compares the suffixes in place.
   The empty suffix (offset == length) compares smallest, which is exactly where '$' sorts. */
int cmpSA(const void *a, const void *b) {
    return strcmp(sortText + *(const int *) a, sortText + *(const int *) b);
}

/* builds suffix array of seq + '$' (length + 1 entries; entry 0 is always length) */
int *buildSuffixArray(char *seq, int length) {
    int *sa = xmalloc(sizeof(int) * (length + 1));
    for (int i = 0; i <= length; i++) {
        sa[i] = i;
    }
    sortText = seq;
    qsort(sa, length + 1, sizeof(int), &cmpSA);
    sortText = NULL;
    return sa;
}

/* builds the suffix array, BWT, F column, and rank tables for fm->seq */
void buildIndex(FM *fm) {
    int n = fm->length;
    int rows = n + 1;
    const unsigned char *s = (const unsigned char *) fm->seq;

    fm->suffixArray = buildSuffixArray(fm->seq, n);

    /* alphabet: '$' is code 0, then every character present, in byte (= sort) order */
    bool present[256] = {false};
    for (int i = 0; i < n; i++) {
        present[s[i]] = true;
    }
    fm->alphabetSize = 1;
    for (int c = 0; c < 256; c++) {
        fm->code[c] = present[c] ? fm->alphabetSize++ : -1;
    }
    int sigma = fm->alphabetSize;

    /* BWT (L column) and F column, straight from the suffix array:
       row r is the rotation starting at sa[r]; L is the character before it, F is the one at it */
    fm->bwt = xmalloc(rows + 1);
    fm->F = xmalloc(rows + 1);
    fm->bwtCodes = xmalloc(rows);
    for (int r = 0; r < rows; r++) {
        int start = fm->suffixArray[r];
        fm->F[r] = (start == n) ? FM_SENTINEL : fm->seq[start];
        if (start == 0) {
            fm->bwt[r] = FM_SENTINEL;
            fm->bwtCodes[r] = 0;
        } else {
            fm->bwt[r] = fm->seq[start - 1];
            fm->bwtCodes[r] = (unsigned char) fm->code[s[start - 1]];
        }
    }
    fm->bwt[rows] = '\0';
    fm->F[rows] = '\0';

    /* C[c]: number of characters in T that sort before c */
    fm->C = xcalloc(sigma, sizeof(int));
    int *counts = xcalloc(sigma, sizeof(int));
    for (int r = 0; r < rows; r++) {
        counts[fm->bwtCodes[r]]++;
    }
    for (int c = 1; c < sigma; c++) {
        fm->C[c] = fm->C[c - 1] + counts[c - 1];
    }
    free(counts);

    /* occurrence checkpoints: occ[k][c] = count of c in bwt[0, k * OCC_SAMPLE) */
    int checkpoints = rows / OCC_SAMPLE + 1;
    fm->occ = xcalloc((size_t) checkpoints * sigma, sizeof(int));
    int *running = xcalloc(sigma, sizeof(int));
    for (int r = 0; r <= rows; r++) {
        if (r % OCC_SAMPLE == 0) {
            memcpy(fm->occ + (size_t) (r / OCC_SAMPLE) * sigma, running, sizeof(int) * sigma);
        }
        if (r < rows) {
            running[fm->bwtCodes[r]]++;
        }
    }
    free(running);
}

/* number of occurrences of code c in bwt[0, row) */
int occCount(FM *fm, int c, int row) {
    int k = row / OCC_SAMPLE;
    int count = fm->occ[(size_t) k * fm->alphabetSize + c];
    for (int r = k * OCC_SAMPLE; r < row; r++) {
        count += (fm->bwtCodes[r] == c);
    }
    return count;
}

/* FM-Index backward search: rows [start, end) whose rotations begin with pattern */
bool backwardSearch(FM *fm, const char *pattern, int m, Interval *interval) {
    int start = 0;
    int end = fm->length + 1;
    for (int i = m - 1; i >= 0 && start < end; i--) {
        int c = fm->code[(unsigned char) pattern[i]];
        if (c < 0) {            /* character never occurs in the reference */
            start = end;
            break;
        }
        start = fm->C[c] + occCount(fm, c, start);
        end = fm->C[c] + occCount(fm, c, end);
    }
    interval->start = start;
    interval->end = end;
    return start < end;
}

/* character at position i of T = seq + '$' */
static char textAt(FM *fm, int i) {
    return (i == fm->length) ? FM_SENTINEL : fm->seq[i];
}

/* write our fm-index to output */
void writeFM(FM *fm, FILE *f) {
    int rows = fm->length + 1;

    fprintf(f, "************\n");
    fprintf(f, "* FM-INDEX *\n");
    fprintf(f, "************\n\n");

    /* name */
    fprintf(f, "> %s\n", fm->name);

    /* sequence length */
    fprintf(f, "length: %d\n\n", fm->length);

    /* sequence, with its end marker */
    fprintf(f, "Sequence\n");
    fprintf(f, "%s%c\n\n", fm->seq, FM_SENTINEL);

    /* suffix array */
    int lineCount = 0;
    fprintf(f, "Suffix Array\n");
    for (int i = 0; i < rows; i++) {
        lineCount++;
        if (i == rows - 1) {
            fprintf(f, "%d\n\n", fm->suffixArray[i]);
        } else if (lineCount == 14) {
            fprintf(f, "%d,\n", fm->suffixArray[i]);
            lineCount = 0;
        } else {
            fprintf(f, "%d, ", fm->suffixArray[i]);
        }
    }

    /* suffixes */
    fprintf(f, "Suffixes\n");
    for (int i = 0; i < rows; i++) {
        fprintf(f, "%d:\t %s%c\n", i, fm->seq + fm->suffixArray[i], FM_SENTINEL);
    }
    fprintf(f, "\n");

    /* bwm, generated one rotation at a time */
    fprintf(f, "Burrows-Wheeler Matrix\n");
    for (int i = 0; i < rows; i++) {
        for (int k = 0; k < rows; k++) {
            fputc(textAt(fm, (fm->suffixArray[i] + k) % rows), f);
        }
        fputc('\n', f);
    }
    fprintf(f, "\n");

    /* bwt */
    fprintf(f, "BWT: ");
    fprintf(f, "%s\n\n", fm->bwt);

    /* F and L */
    fprintf(f, "F: ");
    fprintf(f, "%s\n", fm->F);
    fprintf(f, "L: ");
    fprintf(f, "%s\n", fm->bwt);
}

/* frees the index structures (not name/seq) */
void freeIndex(FM *fm) {
    free(fm->suffixArray);
    free(fm->bwt);
    free(fm->F);
    free(fm->bwtCodes);
    free(fm->C);
    free(fm->occ);
    fm->suffixArray = NULL;
    fm->bwt = fm->F = NULL;
    fm->bwtCodes = NULL;
    fm->C = fm->occ = NULL;
}

/* destroys our FM-Index*/
void destroy(FM *fm) {
    freeIndex(fm);
    free(fm->seq);
    free(fm->name);
    free(fm);
}
