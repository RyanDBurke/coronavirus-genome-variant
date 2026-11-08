#include "misc.h"
#include <getopt.h>
#include <stdarg.h>
#include <unistd.h>

bool quietMode = false;

/* files the presets use */
#define DEFAULT_OUT       "./mappings/mapping.sam"
#define DEFAULT_INDEX_OUT "./fm-output/FMindex.txt"
#define DEMO_REF          "./default/ref-small.fa"
#define DEMO_READS        "./default/reads-small.fa"
#define COVID_REF         "./2019-nCoV.fa"

/* print usage to f */
void usage(FILE *f) {
    fprintf(f,
        "usage: run [options] -r REF.fa -q READS.fa\n"
        "       run [options]                  small demo (default/ref-small.fa)\n"
        "       run [options] covid 1K|10K|1M  SARS-CoV-2 genome with the bundled reads\n"
        "\n"
        "Aligns single-end reads to a reference with an FM-Index seed search and\n"
        "fitting alignment, and writes SAM. FASTA inputs may be gzipped.\n"
        "\n"
        "options:\n"
        "  -r, --ref FILE        reference FASTA (first sequence is used)\n"
        "  -q, --reads FILE      reads FASTA\n"
        "  -o, --out FILE        SAM output, '-' for stdout (default: %s)\n"
        "  -i, --index-out FILE  FM-Index text dump, written for references up to 50 bp\n"
        "                        (default: %s)\n"
        "      --quiet           no progress bar or status messages\n"
        "  -h, --help            show this help\n",
        DEFAULT_OUT, DEFAULT_INDEX_OUT);
}

/* sets the inputs for one of the "covid" presets; returns false for an unknown size */
static bool covidPreset(const char *size, Options *opt) {
    static const struct { const char *name, *reads; double count; } presets[] = {
        { "1K",  "./reads/reads_1K.fa.gz",  1000 },
        { "10K", "./reads/reads_10K.fa.gz", 10000 },
        { "1M",  "./reads/reads_1M.fa.gz",  1000000 },
    };
    for (size_t i = 0; i < sizeof presets / sizeof presets[0]; i++) {
        if (strcasecmp(size, presets[i].name) == 0) {
            opt->ref = COVID_REF;
            opt->reads = presets[i].reads;
            opt->readCount = presets[i].count;
            return true;
        }
    }
    return false;
}

/* parse the command line into opt */
int parseArgs(int argc, char **argv, Options *opt) {
    static const struct option longOptions[] = {
        { "ref",       required_argument, NULL, 'r' },
        { "reads",     required_argument, NULL, 'q' },
        { "out",       required_argument, NULL, 'o' },
        { "index-out", required_argument, NULL, 'i' },
        { "quiet",     no_argument,       NULL, 'Q' },
        { "help",      no_argument,       NULL, 'h' },
        { NULL, 0, NULL, 0 }
    };

    *opt = (Options) { .out = DEFAULT_OUT, .indexOut = DEFAULT_INDEX_OUT };

    int c;
    opterr = 0;  /* we print our own errors */
    while ((c = getopt_long(argc, argv, ":r:q:o:i:h", longOptions, NULL)) != -1) {
        switch (c) {
            case 'r': opt->ref = optarg; break;
            case 'q': opt->reads = optarg; break;
            case 'o': opt->out = optarg; break;
            case 'i': opt->indexOut = optarg; break;
            case 'Q': opt->quiet = true; break;
            case 'h': usage(stdout); return ARGS_HELP;
            case ':':
                fprintf(stderr, "run: option '%s' needs a value\n", argv[optind - 1]);
                return ARGS_ERROR;
            default:
                fprintf(stderr, "run: unknown option '%s'\n", argv[optind - 1]);
                return ARGS_ERROR;
        }
    }

    int positional = argc - optind;
    bool explicitInputs = opt->ref != NULL || opt->reads != NULL;

    if (positional == 0 && !explicitInputs) {
        /* [run] small demo */
        opt->ref = DEMO_REF;
        opt->reads = DEMO_READS;
        opt->readCount = 1;
    } else if (positional == 2 && !explicitInputs && strcasecmp(argv[optind], "covid") == 0) {
        /* [run covid 1K|10K|1M] */
        if (!covidPreset(argv[optind + 1], opt)) {
            fprintf(stderr, "run: unknown preset size '%s' (expected 1K, 10K or 1M)\n", argv[optind + 1]);
            return ARGS_ERROR;
        }
    } else if (positional == 0) {
        /* [run -r REF -q READS] */
        if (opt->ref == NULL || opt->reads == NULL) {
            fprintf(stderr, "run: --ref and --reads must be given together\n");
            return ARGS_ERROR;
        }
    } else {
        fprintf(stderr, "run: unexpected argument '%s'\n", argv[optind]);
        return ARGS_ERROR;
    }

    /* SAM on stdout must not be mixed with status output */
    if (strcmp(opt->out, "-") == 0) {
        opt->quiet = true;
    }
    return ARGS_OK;
}

/* status messages */
void info(const char *fmt, ...) {
    if (quietMode) {
        return;
    }
    va_list args;
    va_start(args, fmt);
    vfprintf(stderr, fmt, args);
    va_end(args);
}

const char *highlight(void) {
    return isatty(STDERR_FILENO) ? "\033[1;31m" : "";
}

const char *highlightEnd(void) {
    return isatty(STDERR_FILENO) ? "\033[0m" : "";
}

/* progress bar: redraws only when the whole-number percentage changes */
void pb(bool makeProgressBar, double percentageDone) {
    static int lastPercent = -1;
    if (!makeProgressBar || quietMode) {
        return;
    }

    /* clamp to [0, 100] */
    int percent = (int) (percentageDone * 100);
    percent = percent < 0 ? 0 : (percent > 100 ? 100 : percent);
    if (percent == lastPercent) {
        return;
    }
    lastPercent = percent;

    /* 12-character bar with at least one '#' */
    int filled = 1 + percent * 11 / 100;
    fprintf(stderr, "\rProgress: [%.*s%*s] %d%%", filled, "############", 12 - filled, "", percent);
    fflush(stderr);
}

/* allocation wrappers: report and exit instead of returning NULL */
static void outOfMemory(size_t bytes) {
    fprintf(stderr, "out of memory allocating %zu bytes\n", bytes);
    exit(1);
}

void *xmalloc(size_t size) {
    void *p = malloc(size);
    if (p == NULL && size != 0) {
        outOfMemory(size);
    }
    return p;
}

void *xcalloc(size_t count, size_t size) {
    void *p = calloc(count, size);
    if (p == NULL && count != 0 && size != 0) {
        outOfMemory(count * size);
    }
    return p;
}

void *xrealloc(void *ptr, size_t size) {
    void *p = realloc(ptr, size);
    if (p == NULL && size != 0) {
        outOfMemory(size);
    }
    return p;
}
