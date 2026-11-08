#include "align.h"
#include "fmIndex.h"
#include "misc.h"
#include <errno.h>
#include <limits.h>
#include <sys/stat.h>
KSEQ_INIT(gzFile, gzread)

/* SAM flags we emit */
#define FLAG_UNMAPPED  0x4
#define FLAG_SECONDARY 0x100

/* longest read the fixed-size OPT-matrix can hold */
#define MAX_READ_LENGTH (MAXCOL - 2 * SLICE_PAD - 1)

/* write a record for a read with no alignment */
static void writeUnmapped(FILE *f, kseq_t *seq) {
    SAM s = {
        .QNAME = seq->name.s, .FLAG = FLAG_UNMAPPED, .RNAME = "*", .POS = 0, .MAPQ = 0,
        .CIGAR = "*", .RNEXT = "*", .PNEXT = 0, .TLEN = 0, .SEQ = seq->seq.s, .QUAL = "*"
    };
    writeSAM(&s, f);
}

/* size of a file in bytes, or 0 if unknown */
static double fileSize(const char *path) {
    struct stat st;
    return (stat(path, &st) == 0 && st.st_size > 0) ? (double) st.st_size : 0;
}

int align(FM *fm, const char *reads, const char *output, AlignStats *stats) {

    *stats = (AlignStats) {0};

    /* fasta stuff */
    if (reads == NULL) {
        return fail(EXIT_USAGE, "no reads file given");
    }
    gzFile fp = gzopen(reads, "rb");
    if (fp == NULL) {
        return fail(EXIT_INPUT, "can't open reads file '%s': %s", reads, strerror(errno));
    }

    /* progress is the share of the (possibly compressed) file read so far, so it works for any input */
    double readsBytes = fileSize(reads);

    /* open our output file (.sam), or write to stdout for "-"; binary keeps "\n" line endings on Windows */
    bool toStdout = strcmp(output, "-") == 0;
    FILE *f = toStdout ? stdout : fopen(output, "wb");
    if (f == NULL) {
        int code = fail(EXIT_OUTPUT, "can't write output file '%s': %s", output, strerror(errno));
        gzclose(fp);
        return code;
    }

    kseq_t *seq = kseq_init(fp);
    int gap = -5; /* gap penalty */
    bool warnedLongRead = false;

    /* .sam header */
    fprintf(f, "@HD\tVN:1.0\tSO:unsorted\n");
    fprintf(f, "@SQ\tSN:%s\tLN:%d\n", fm->name, fm->length);

    /* best-scoring alignments of the current read, across all of its seeds */
    Alignments hits = {0};
    info("> Aligning reads from %s\n", reads);

    /* parse .fa file containing our reads */
    int readStatus;
    while ((readStatus = kseq_read(seq)) >= 0) {

        stats->reads++;
        if (readsBytes > 0) {
            progress((double) gzoffset(fp) / readsBytes);
        }

        int readLength = (int) seq->seq.l;
        int bestScore = INT_MIN;
        clearAlignments(&hits);

        if (readLength > MAX_READ_LENGTH) {
            if (!warnedLongRead) {
                warn("reads longer than %d bp are reported as unmapped", MAX_READ_LENGTH);
                warnedLongRead = true;
            }
        } else if (readLength > 0) {

            /* for each (read-length / 5) seed (20bp seeds for our 100bp reads) */
            int skip = seedSkip(readLength);
            char *seed = xmalloc(skip + 1);
            for (int seedStart = 0; seedStart < readLength; seedStart += skip) {
                int seedEnd = min(readLength, seedStart + skip);

                /* find where this seed matches the reference exactly */
                Interval interval;
                int matchLength = 0;
                substring(seed, seq->seq.s, seedStart, seedEnd);
                getInterval(&interval, &matchLength, fm, seed);

                /* no exact match (e.g. a variant inside this seed): try the next seed */
                if (matchLength == 0) {
                    continue;
                }

                /* candidate read start positions in the reference */
                int *refPos = xmalloc(sizeof(int) * (interval.end - interval.start));
                int refPosLength = referencePos(refPos, &interval, matchLength, fm, seedEnd);

                /* fitting alignment at each candidate, keeping the best across all seeds */
                alignment(&hits, seq->seq.s, readLength, fm->seq, fm->length, refPos, refPosLength, gap, &bestScore);

                free(refPos);
            }
            free(seed);
        }

        /* write to .sam file: first hit is primary, any equally good ones are secondary */
        if (hits.count == 0) {
            writeUnmapped(f, seq);
        } else {
            stats->mapped++;
        }
        for (int a = 0; a < hits.count; a++) {
            SAM s = {
                .QNAME = seq->name.s,
                .FLAG = (a == 0) ? 0 : FLAG_SECONDARY,
                .RNAME = fm->name,
                .POS = hits.items[a].pos,
                .MAPQ = 255,                 /* not computed */
                .CIGAR = hits.items[a].cigar,
                .RNEXT = "*",
                .PNEXT = 0,
                .TLEN = 0,                   /* single-end reads */
                .SEQ = seq->seq.s,
                .QUAL = "*"
            };
            writeSAM(&s, f);
        }
    }

    clearAlignments(&hits);
    free(hits.items);
    kseq_destroy(seq);
    gzclose(fp);

    /* a failed write (e.g. a full disk) only shows up here */
    bool writeFailed = ferror(f) != 0;
    writeFailed |= toStdout ? fflush(f) != 0 : fclose(f) != 0;
    if (writeFailed) {
        return fail(EXIT_OUTPUT, "error writing output file '%s'", toStdout ? "stdout" : output);
    }

    /* -1 is a clean end of file; less means a truncated or corrupt file */
    if (readStatus < -1) {
        return fail(EXIT_INPUT, "reads file '%s' is truncated or corrupt", reads);
    }
    if (stats->reads == 0) {
        return fail(EXIT_INPUT, "no reads found in '%s' (expected FASTA or FASTQ)", reads);
    }

    progress(1.0);
    info("\n");  /* end the progress bar line */
    info("> %d of %d reads mapped; SAM alignments written to %s%s%s\n",
         stats->mapped, stats->reads, highlight(), output, highlightEnd());
    return EXIT_OK;
}

/**************************/
/* AUX FUNTIONS FOR ALIGN */
/**************************/

/* finds the suffix-array interval [start, end) whose suffixes begin with seed (FM-Index backward search) */
void getInterval(Interval *interval, int *matchLength, FM *fm, char* seed) {
    int seedLength = strlen(seed);

    if (seedLength > 0 && backwardSearch(fm, seed, seedLength, interval)) {
        *matchLength = seedLength;
    } else {
        interval->start = -1;
        interval->end = -1;
        *matchLength = 0;
    }
}

/* returns reference positions of current seed */
int referencePos(int *refPos, Interval *interval, int matchLength, FM *fm, int seedEnd) {

    /* interval: [start, end) */
    int start = interval->start;
    int end = interval->end;

    /* the read starts (seedEnd - matchLength) bases before where the seed matched */
    int refPosLength = 0;
    for (int i = start; i < end; i++) {
        int pos = (fm->suffixArray[i]) - (seedEnd - matchLength);
        refPos[refPosLength] = pos;
        refPosLength++;
    }

    /* we're returning the length of refPos int-array */
    return refPosLength;
}

/* empties an alignment list, freeing its CIGAR strings */
void clearAlignments(Alignments *hits) {
    for (int i = 0; i < hits->count; i++) {
        free(hits->items[i].cigar);
    }
    hits->count = 0;
}

/* appends an alignment, growing the list as needed; takes ownership of cigar */
static void addAlignment(Alignments *hits, int score, int pos, char *cigar) {
    if (hits->count == hits->capacity) {
        hits->capacity = hits->capacity ? hits->capacity * 2 : 4;
        hits->items = xrealloc(hits->items, sizeof(Alignment) * hits->capacity);
    }
    hits->items[hits->count].score = score;
    hits->items[hits->count].pos = pos;
    hits->items[hits->count].cigar = cigar;
    hits->count++;
}

/* fitting alignment of read at each candidate position; keeps the best-scoring hits */
void alignment(Alignments *hits, char *read, int readLength, char *ref, int refLength, int *refPos, int refPosLength, int gap, int *bestScore) {

    /* X-axis is our coronavirus reference genome slice */
    /* Y-axis is our read we wish to align to our reference */
    int OPT[MAXROW][MAXCOL];
    char *x = xmalloc(readLength + 2 * SLICE_PAD + 1);

    for (int p = 0; p < refPosLength; p++) {

        /* slice out the reference with SLICE_PAD bases on each side, clipped to the genome */
        int sliceStart = max(refPos[p] - SLICE_PAD, 0);
        int sliceEnd = min(refPos[p] + readLength + SLICE_PAD, refLength);
        if (sliceStart >= sliceEnd) {
            continue;
        }
        substring(x, ref, sliceStart, sliceEnd);
        int n = sliceEnd - sliceStart;

        /* only keep alignments with >= best score */
        int score = editDistance(OPT, x, read, n, readLength, gap);
        if (score < *bestScore) {
            continue;
        }
        if (score > *bestScore) {
            clearAlignments(hits);
            *bestScore = score;
        }

        /* backtrace OPT-matrix to find alignment and its CIGAR string */
        int offset = -1;
        char *cigar = buildCigar(OPT, n, readLength, gap, x, read, &offset);

        /* offset is where the read starts inside the slice; SAM POS is 1-based */
        int pos = sliceStart + offset + 1;

        /* different seeds of one read usually find the same locus; keep it once */
        bool duplicate = false;
        for (int a = 0; a < hits->count; a++) {
            if (hits->items[a].pos == pos) {
                duplicate = true;
                break;
            }
        }
        if (duplicate) {
            free(cigar);
        } else {
            addAlignment(hits, score, pos, cigar);
        }
    }

    free(x);
}

/* builds OPT-matrix for a fitting alignment of y inside x, and returns the best score
   (max over the last column: reference bases before or after the read are free) */
int editDistance(int OPT[MAXROW][MAXCOL], char *x, char *y, int n, int m, int gap) {

    /* (0, 0) is always 0 */
    OPT[0][0] = 0;

    /* cost-free gaps before the read starts on the reference */
    for (int i = 1; i <= n; i++) {
        OPT[i][0] = 0;
    }

    /* add our initial gap penalties to the first row of each column */
    for (int j = 1; j <= m; j++) {
        OPT[0][j] = j * gap;
    }

    /* fill matrix */
    for (int i = 1; i <= n; i++) {
        for (int j = 1; j <= m; j++) {
            OPT[i][j] = maxAlign(
                    (score(x[i - 1], y[j - 1]) + OPT[i - 1][j - 1]),
                    (gap + OPT[i - 1][j]),
                    (gap + OPT[i][j - 1])
                );
        }
    }

    /* the read may end anywhere on the slice */
    int best = OPT[0][m];
    for (int i = 1; i <= n; i++) {
        best = max(best, OPT[i][m]);
    }
    return best;
}

/* returns CIGAR-string */
char *buildCigar(int OPT[MAXROW][MAXCOL], int n, int m, int gap, char *x, char *y, int *offset) {

    /* the read ends at the best-scoring row of the last column */
    int row = 0;
    int maxScore = INT_MIN;
    for (int i = 0; i <= n; i++) {
        if (OPT[i][m] >= maxScore) {
            maxScore = OPT[i][m];
            row = i;
        }
    }
    n = row;

    /* backtrace from the end of the read to its start, collecting one op per step */
    char *ops = xmalloc(n + m + 1);
    int length = 0;
    while (m > 0) {
        if (n > 0 && OPT[n][m] == score(x[n - 1], y[m - 1]) + OPT[n - 1][m - 1]) {
            ops[length++] = 'M';
            n--;
            m--;
        } else if (OPT[n][m] == gap + OPT[n][m - 1]) {
            ops[length++] = 'I';
            m--;
        } else if (n > 0 && OPT[n][m] == gap + OPT[n - 1][m]) {
            ops[length++] = 'D';
            n--;
        } else {
            exit(fail(EXIT_GENERIC, "error in buildCigar(): inconsistent OPT-matrix"));
        }
    }

    /* we stopped at the slice position where the read starts */
    *offset = n;

    /* ops were collected end-to-start; reverse them into read order */
    for (int i = 0, j = length - 1; i < j; i++, j--) {
        char t = ops[i];
        ops[i] = ops[j];
        ops[j] = t;
    }
    ops[length] = '\0';

    /* run-length encode into a CIGAR string (i.e MMMDDMM -> 3M2D2M) */
    char *cigar = formatCigar(ops);
    free(ops);
    return cigar;
}

/* return score between two characters */
int score(char a, char b) {
    if (a == 'N' || b == 'N') { /* if you see 'N', assume it matches */
        return 0;
    }
    return (a == b) ? 1 : -1;
}

/* return substring of string from [start, end) */
void substring(char* result, char* string, int start, int end) {
    int n = end - start;
    memcpy(result, string + start, n);
    result[n] = '\0';
}

/* return correctly formatted CIGAR-string (i.e MMMDDMM -> 3M2D2M) */
char *formatCigar(const char *ops) {

    /* a run of length r takes at most r + 1 characters, so 2 * length is enough */
    size_t length = strlen(ops);
    char *result = xmalloc(2 * length + 1);
    size_t written = 0;

    for (size_t i = 0; i < length;) {
        size_t j = i;
        while (j < length && ops[j] == ops[i]) {
            j++;
        }
        char run[32];
        int runLength = snprintf(run, sizeof run, "%zu%c", j - i, ops[i]);
        memcpy(result + written, run, runLength);
        written += runLength;
        i = j;
    }
    result[written] = '\0';

    return result;
}

/* seed skip (and seed length): a fifth of the read, at least 1 */
int seedSkip(int L) {
    return max(L / 5, 1);
}

/* min */
int min(int a, int b) {
    return (a > b) ? b : a;
}

/* max */
int max(int a, int b) {
    return (a > b) ? a : b;
}

/* max function needed to fill OPT-matrix */
int maxAlign(int a, int b, int c) {
    return max(a, max(b, c));
}

/* write to our .SAM file */
void writeSAM(SAM *sam, FILE *f) {
    fprintf(f, "%s\t%d\t%s\t%d\t%d\t%s\t%s\t%d\t%d\t%s\t%s\n",
            sam->QNAME, sam->FLAG, sam->RNAME, sam->POS, sam->MAPQ, sam->CIGAR,
            sam->RNEXT, sam->PNEXT, sam->TLEN, sam->SEQ, sam->QUAL);
}
