/* Unit tests for the alignment and indexing helpers.
 * Built with AddressSanitizer + UBSan by `make test`. */

#include "align.h"
#include "fmIndex.h"

static int failures = 0;
static int checks = 0;

#define CHECK(cond, ...) do {                                          \
        checks++;                                                       \
        if (!(cond)) {                                                  \
            failures++;                                                 \
            printf("FAIL %s:%d: ", __FILE__, __LINE__);                 \
            printf(__VA_ARGS__);                                        \
            printf("\n");                                               \
        }                                                               \
    } while (0)

/* reference context with a unique 20bp core; padding mimics the 5bp slice padding */
#define CORE "CATGCTTAGTGCACTCACGC"
#define PADDED "AAAAA" CORE "AAAAA"

/* runs the fitting alignment of read y against slice x and returns its CIGAR */
static char *alignPair(const char *x, const char *y, int *offset) {
    static int OPT[MAXROW][MAXCOL];
    int n = strlen(x), m = strlen(y);
    editDistance(OPT, (char *)x, (char *)y, n, m, -5);
    *offset = -1;
    return buildCigar(OPT, n, m, -5, (char *)x, (char *)y, offset);
}

/* number of read bases a CIGAR consumes (M + I) */
static int queryLength(const char *cigar) {
    int total = 0, n = 0;
    for (const char *c = cigar; *c; c++) {
        if (*c >= '0' && *c <= '9') {
            n = n * 10 + (*c - '0');
        } else {
            if (*c == 'M' || *c == 'I') total += n;
            n = 0;
        }
    }
    return total;
}

static void expectCigar(const char *label, const char *x, const char *y, const char *cigar, int offset) {
    int got_offset;
    char *got = alignPair(x, y, &got_offset);
    CHECK(strcmp(got, cigar) == 0, "%s: CIGAR %s, expected %s", label, got, cigar);
    CHECK(got_offset == offset, "%s: offset %d, expected %d", label, got_offset, offset);
    free(got);
}

static void test_cigars(void) {
    /* exact match: 20M starting 5 bases into the slice */
    expectCigar("exact", PADDED, CORE, "20M", 5);

    /* a mismatch is still M */
    expectCigar("snp", PADDED, "CATGCATAGTGCACTCACGC", "20M", 5);

    /* read has an extra G between CORE[11] and CORE[12] */
    expectCigar("insertion", PADDED, "CATGCTTAGTGCGACTCACGC", "12M1I8M", 5);

    /* read is missing CORE[10]; asymmetric so a reversed CIGAR is caught */
    expectCigar("deletion", PADDED, "CATGCTTAGTCACTCACGC", "10M1D9M", 5);

    /* runs of 10+ need multi-digit counts */
    expectCigar("long run", "AAAAA" CORE CORE "AAAAA", CORE CORE, "40M", 5);
}

static void test_poor_alignment_stays_in_bounds(void) {
    /* every score in the last column is below -1; used to read OPT[-1] */
    int offset;
    char *cigar = alignPair("AAAA", "CCCCCCCC", &offset);
    CHECK(queryLength(cigar) == 8, "poor alignment: CIGAR %s should consume 8 read bases", cigar);
    CHECK(offset >= 0 && offset <= 4, "poor alignment: offset %d out of range", offset);
    free(cigar);
}

static void test_interval_reaching_end_of_suffix_array(void) {
    /* rows of ACGTTT$ sorted: $ ACGTTT$ CGTTT$ GTTT$ T$ TT$ TTT$
       seed "TT" matches the last two rows [5, 7) */
    char seq[] = "ACGTTT";
    FM fm = {0};
    fm.seq = seq;
    fm.length = strlen(seq);
    buildIndex(&fm);

    CHECK(strcmp(fm.bwt, "T$ACTTG") == 0, "BWT of ACGTTT: %s, expected T$ACTTG", fm.bwt);
    CHECK(strcmp(fm.F, "$ACGTTT") == 0, "F of ACGTTT: %s, expected $ACGTTT", fm.F);

    Interval interval;
    int matchLength = 0;
    getInterval(&interval, &matchLength, &fm, "TT");
    CHECK(matchLength == 2, "interval at end: matchLength %d, expected 2", matchLength);
    CHECK(interval.start == 5 && interval.end == 7,
          "interval at end: [%d, %d), expected [5, 7)", interval.start, interval.end);

    getInterval(&interval, &matchLength, &fm, "CG");
    CHECK(interval.start == 2 && interval.end == 3,
          "interval in middle: [%d, %d), expected [2, 3)", interval.start, interval.end);

    getInterval(&interval, &matchLength, &fm, "GG");
    CHECK(matchLength == 0, "no match: matchLength %d, expected 0", matchLength);

    getInterval(&interval, &matchLength, &fm, "TN");
    CHECK(matchLength == 0, "character absent from reference: matchLength %d, expected 0", matchLength);

    freeIndex(&fm);
}

/* backward search must find exactly the positions a brute-force scan finds */
static void test_backward_search_matches_brute_force(void) {
    /* pseudo-random sequence over ACGT plus a few N, with repeats so patterns hit many times */
    enum { LENGTH = 3000 };
    char *seq = malloc(LENGTH + 1);
    unsigned int state = 2020;
    for (int i = 0; i < LENGTH; i++) {
        state = state * 1103515245u + 12345u;
        seq[i] = (i % 997 == 0) ? 'N' : "ACGT"[(state >> 16) % 4];
    }
    memcpy(seq + 1500, seq + 100, 200);  /* a long repeat */
    seq[LENGTH] = '\0';

    FM fm = {0};
    fm.seq = seq;
    fm.length = LENGTH;
    buildIndex(&fm);

    int mismatches = 0;
    int total = 0;
    for (int trial = 0; trial < 400; trial++) {
        /* patterns cut from the sequence (always present) and random ones (often absent) */
        int m = 1 + trial % 25;
        char pattern[32];
        if (trial % 2 == 0) {
            int at = (trial * 7919) % (LENGTH - m);
            memcpy(pattern, seq + at, m);
        } else {
            for (int k = 0; k < m; k++) {
                state = state * 1103515245u + 12345u;
                pattern[k] = "ACGT"[(state >> 16) % 4];
            }
        }
        pattern[m] = '\0';

        /* brute force occurrence count */
        int expected = 0;
        for (int i = 0; i + m <= LENGTH; i++) {
            expected += (memcmp(seq + i, pattern, m) == 0);
        }

        Interval interval;
        bool found = backwardSearch(&fm, pattern, m, &interval);
        int got = found ? interval.end - interval.start : 0;

        /* every row in the interval must really start with the pattern */
        bool rowsOk = true;
        for (int r = interval.start; found && r < interval.end; r++) {
            int at = fm.suffixArray[r];
            rowsOk &= (at + m <= LENGTH && memcmp(seq + at, pattern, m) == 0);
        }
        if (got != expected || !rowsOk) {
            mismatches++;
        }
        total++;
    }
    CHECK(mismatches == 0, "backward search disagreed with brute force on %d of %d patterns", mismatches, total);

    freeIndex(&fm);
    free(seq);
}

static void test_seed_skip(void) {
    CHECK(seedSkip(100) == 20, "seedSkip(100) = %d, expected 20", seedSkip(100));
    CHECK(seedSkip(3) >= 1, "seedSkip(3) = %d; must be >= 1 or the seed loop never advances", seedSkip(3));
}

static void test_score(void) {
    CHECK(score('A', 'A') == 1, "match should score 1");
    CHECK(score('A', 'C') == -1, "mismatch should score -1");
    CHECK(score('N', 'C') == 0, "N should score 0");
}

int main(void) {
    /* unbuffered, so FAIL lines survive a sanitizer abort */
    setvbuf(stdout, NULL, _IONBF, 0);

    test_cigars();
    test_poor_alignment_stays_in_bounds();
    test_interval_reaching_end_of_suffix_array();
    test_backward_search_matches_brute_force();
    test_seed_skip();
    test_score();

    printf("unit tests: %d checks, %d failures\n", checks, failures);
    return failures == 0 ? 0 : 1;
}
