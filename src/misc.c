#include "misc.h"
#include <getopt.h>
#include <stdarg.h>
#include <unistd.h>
#ifdef _WIN32
#include <windows.h>
#endif

bool quietMode = false;
bool jsonMode = false;

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
        "      --progress-json   machine-readable progress: one JSON object per line on stderr\n"
        "  -h, --help            show this help\n"
        "\n"
        "exit codes: 0 ok, 1 internal error, 2 invalid usage, 3 input file problem,\n"
        "            4 output file problem, 5 out of memory\n",
        DEFAULT_OUT, DEFAULT_INDEX_OUT);
}

/* sets the inputs for one of the "covid" presets; returns false for an unknown size */
static bool covidPreset(const char *size, Options *opt) {
    static const struct { const char *name, *reads; } presets[] = {
        { "1K",  "./reads/reads_1K.fa.gz" },
        { "10K", "./reads/reads_10K.fa.gz" },
        { "1M",  "./reads/reads_1M.fa.gz" },
    };
    for (size_t i = 0; i < sizeof presets / sizeof presets[0]; i++) {
        if (strcasecmp(size, presets[i].name) == 0) {
            opt->ref = COVID_REF;
            opt->reads = presets[i].reads;
            return true;
        }
    }
    return false;
}

/* parse the command line into opt */
int parseArgs(int argc, char **argv, Options *opt) {
    static const struct option longOptions[] = {
        { "ref",           required_argument, NULL, 'r' },
        { "reads",         required_argument, NULL, 'q' },
        { "out",           required_argument, NULL, 'o' },
        { "index-out",     required_argument, NULL, 'i' },
        { "quiet",         no_argument,       NULL, 'Q' },
        { "progress-json", no_argument,       NULL, 'J' },
        { "help",          no_argument,       NULL, 'h' },
        { NULL, 0, NULL, 0 }
    };

    *opt = (Options) { .out = DEFAULT_OUT, .indexOut = DEFAULT_INDEX_OUT };

    /* JSON mode affects how argument errors are reported, so find it first */
    for (int i = 1; i < argc; i++) {
        if (strcmp(argv[i], "--progress-json") == 0) {
            jsonMode = true;
        } else if (strcmp(argv[i], "--") == 0) {
            break;
        }
    }

    int c;
    opterr = 0;  /* we print our own errors */
    while ((c = getopt_long(argc, argv, ":r:q:o:i:h", longOptions, NULL)) != -1) {
        switch (c) {
            case 'r': opt->ref = optarg; break;
            case 'q': opt->reads = optarg; break;
            case 'o': opt->out = optarg; break;
            case 'i': opt->indexOut = optarg; break;
            case 'Q': opt->quiet = true; break;
            case 'J': opt->progressJson = true; break;
            case 'h': usage(stdout); return ARGS_HELP;
            case ':':
                fail(EXIT_USAGE, "option '%s' needs a value", argv[optind - 1]);
                return ARGS_ERROR;
            default:
                fail(EXIT_USAGE, "unknown option '%s'", argv[optind - 1]);
                return ARGS_ERROR;
        }
    }

    int positional = argc - optind;
    bool explicitInputs = opt->ref != NULL || opt->reads != NULL;

    if (positional == 0 && !explicitInputs) {
        /* [run] small demo */
        opt->ref = DEMO_REF;
        opt->reads = DEMO_READS;
    } else if (positional == 2 && !explicitInputs && strcasecmp(argv[optind], "covid") == 0) {
        /* [run covid 1K|10K|1M] */
        if (!covidPreset(argv[optind + 1], opt)) {
            fail(EXIT_USAGE, "unknown preset size '%s' (expected 1K, 10K or 1M)", argv[optind + 1]);
            return ARGS_ERROR;
        }
    } else if (positional == 0) {
        /* [run -r REF -q READS] */
        if (opt->ref == NULL || opt->reads == NULL) {
            fail(EXIT_USAGE, "--ref and --reads must be given together");
            return ARGS_ERROR;
        }
    } else {
        fail(EXIT_USAGE, "unexpected argument '%s'", argv[optind]);
        return ARGS_ERROR;
    }

    /* SAM on stdout must not be mixed with status output; JSON events replace it */
    if (strcmp(opt->out, "-") == 0 || opt->progressJson) {
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

/* writes s as a JSON string literal */
static void jsonString(FILE *f, const char *s) {
    fputc('"', f);
    for (const unsigned char *p = (const unsigned char *) s; *p; p++) {
        switch (*p) {
            case '"':  fputs("\\\"", f); break;
            case '\\': fputs("\\\\", f); break;
            case '\n': fputs("\\n", f); break;
            case '\r': fputs("\\r", f); break;
            case '\t': fputs("\\t", f); break;
            default:
                if (*p < 0x20) {
                    fprintf(f, "\\u%04x", *p);
                } else {
                    fputc(*p, f);  /* UTF-8 passes through unchanged */
                }
        }
    }
    fputc('"', f);
}

void event(const char *name, const char *fields, ...) {
    if (!jsonMode) {
        return;
    }
    fprintf(stderr, "{\"event\":\"%s\"", name);
    if (fields[0] != '\0') {
        fputc(',', stderr);
        va_list args;
        va_start(args, fields);
        vfprintf(stderr, fields, args);
        va_end(args);
    }
    fputs("}\n", stderr);
    fflush(stderr);
}

/* human message on stderr, or a JSON event carrying it as "message" */
static void report(const char *eventName, const char *codeName, const char *prefix,
                   const char *fmt, va_list args) {
    char message[1024];
    vsnprintf(message, sizeof message, fmt, args);
    if (jsonMode) {
        fprintf(stderr, "{\"event\":\"%s\",", eventName);
        if (codeName != NULL) {
            fprintf(stderr, "\"code\":\"%s\",", codeName);
        }
        fputs("\"message\":", stderr);
        jsonString(stderr, message);
        fputs("}\n", stderr);
        fflush(stderr);
    } else {
        fprintf(stderr, "%s%s\n", prefix, message);
    }
}

static const char *exitCodeName(int code) {
    switch (code) {
        case EXIT_USAGE:  return "usage";
        case EXIT_INPUT:  return "input";
        case EXIT_OUTPUT: return "output";
        case EXIT_NOMEM:  return "memory";
        default:          return "internal";
    }
}

int fail(int code, const char *fmt, ...) {
    va_list args;
    va_start(args, fmt);
    report("error", exitCodeName(code), "run: ", fmt, args);
    va_end(args);
    return code;
}

void warn(const char *fmt, ...) {
    va_list args;
    va_start(args, fmt);
    report("warning", NULL, "run: warning: ", fmt, args);
    va_end(args);
}

/* true when stderr is a terminal that understands ANSI colour codes */
static bool colourTerminal(void) {
    static int cached = -1;
    if (cached < 0) {
        cached = isatty(STDERR_FILENO) ? 1 : 0;
#ifdef _WIN32
        /* Windows consoles need virtual-terminal processing switched on first */
        if (cached) {
            HANDLE h = GetStdHandle(STD_ERROR_HANDLE);
            DWORD mode = 0;
            cached = GetConsoleMode(h, &mode)
                  && SetConsoleMode(h, mode | ENABLE_VIRTUAL_TERMINAL_PROCESSING);
        }
#endif
    }
    return cached == 1;
}

const char *highlight(void) {
    return colourTerminal() ? "\033[1;31m" : "";
}

const char *highlightEnd(void) {
    return colourTerminal() ? "\033[0m" : "";
}

/* progress: reports only when the whole-number percentage changes */
void progress(double fraction) {
    static int lastPercent = -1;

    /* clamp to [0, 100] */
    int percent = (int) (fraction * 100);
    percent = percent < 0 ? 0 : (percent > 100 ? 100 : percent);
    if (percent == lastPercent) {
        return;
    }
    lastPercent = percent;

    if (jsonMode) {
        event("progress", "\"percent\":%d", percent);
        return;
    }
    if (quietMode) {
        return;
    }

    /* 12-character bar with at least one '#' */
    int filled = 1 + percent * 11 / 100;
    fprintf(stderr, "\rProgress: [%.*s%*s] %d%%", filled, "############", 12 - filled, "", percent);
    fflush(stderr);
}

/* allocation wrappers: report and exit instead of returning NULL */
static void outOfMemory(size_t bytes) {
    exit(fail(EXIT_NOMEM, "out of memory allocating %zu bytes", bytes));
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
