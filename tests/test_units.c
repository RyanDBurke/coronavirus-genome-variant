/* Unit tests for the alignment and indexing helpers.
 * Built with AddressSanitizer + UBSan by `make test`. */

#include "../lib/align.h"
#include "../lib/fmIndex.h"

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
    /* suffixes of ACGTTT sorted: ACGTTT CGTTT GTTT T TT TTT
       seed "TT" matches the last two entries [4, 6) */
    char seq[] = "ACGTTT";
    FM fm = {0};
    fm.seq = seq;
    fm.length = strlen(seq);
    fm.suffixArray = buildSuffixArray(seq, fm.length);

    Interval interval;
    int matchLength = 0;
    getInterval(&interval, &matchLength, &fm, "TT");
    CHECK(matchLength == 2, "interval at end: matchLength %d, expected 2", matchLength);
    CHECK(interval.start == 4 && interval.end == 6,
          "interval at end: [%d, %d), expected [4, 6)", interval.start, interval.end);

    getInterval(&interval, &matchLength, &fm, "CG");
    CHECK(interval.start == 1 && interval.end == 2,
          "interval in middle: [%d, %d), expected [1, 2)", interval.start, interval.end);

    getInterval(&interval, &matchLength, &fm, "GG");
    CHECK(matchLength == 0, "no match: matchLength %d, expected 0", matchLength);

    free(fm.suffixArray);
}

static void test_seed_skip(void) {
    CHECK(seedSkip(100) == 20, "seedSkip(100) = %d, expected 20", seedSkip(100));
    CHECK(seedSkip(3) >= 1, "seedSkip(3) = %d; must be >= 1 or the seed loop never advances", seedSkip(3));
}

static void test_score(void) {
    CHECK(score('A', 'A', -5) == 1, "match should score 1");
    CHECK(score('A', 'C', -5) == -1, "mismatch should score -1");
    CHECK(score('N', 'C', -5) == 0, "N should score 0");
}

int main(void) {
    /* unbuffered, so FAIL lines survive a sanitizer abort */
    setvbuf(stdout, NULL, _IONBF, 0);

    test_cigars();
    test_poor_alignment_stays_in_bounds();
    test_interval_reaching_end_of_suffix_array();
    test_seed_skip();
    test_score();

    printf("unit tests: %d checks, %d failures\n", checks, failures);
    return failures == 0 ? 0 : 1;
}
