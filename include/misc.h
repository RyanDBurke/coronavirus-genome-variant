/****************************************************/
/*                     MISC                         */
/* command-line handling, messages, and allocation */
/****************************************************/

#ifndef MISC_H
#define MISC_H

#include "std.h"

/* printf-style format checking; MinGW needs gnu_printf to accept C99 formats such as %zu (UCRT supports them) */
#ifdef __MINGW32__
#define PRINTF_LIKE(fmt, args) __attribute__((format(gnu_printf, fmt, args)))
#else
#define PRINTF_LIKE(fmt, args) __attribute__((format(printf, fmt, args)))
#endif

/* everything main() needs to know from the command line */
typedef struct options {
    const char  *ref;           // reference FASTA
    const char  *reads;         // reads FASTA
    const char  *out;           // SAM output path, "-" for stdout
    const char  *indexOut;      // FM-Index text dump (references up to 50 bp)
    bool        quiet;          // no progress bar or status messages
    bool        progressJson;   // machine-readable events on stderr instead of human output
} Options;

/* parseArgs results */
#define ARGS_OK    0
#define ARGS_ERROR 1
#define ARGS_HELP  2

/* process exit codes (documented in the README; the desktop app maps them to messages) */
#define EXIT_OK      0
#define EXIT_GENERIC 1   /* internal error */
#define EXIT_USAGE   2   /* invalid command line */
#define EXIT_INPUT   3   /* reference or reads file missing, unreadable or empty */
#define EXIT_OUTPUT  4   /* SAM or FM-Index file could not be written */
#define EXIT_NOMEM   5   /* out of memory */

/* parse the command line into opt; prints its own error messages
    * returns ARGS_OK, ARGS_ERROR (invalid usage) or ARGS_HELP (--help was printed)
 */
int parseArgs(int argc, char **argv, Options *opt);

/* print usage to f */
void usage(FILE *f);

/* status messages go to stderr (stdout may be carrying SAM), and are silenced by --quiet */
extern bool quietMode;
void info(const char *fmt, ...) PRINTF_LIKE(1, 2);

/* --progress-json: every stderr line is one JSON object */
extern bool jsonMode;

/* print one JSON event line: fields is the object body after "event", e.g. "\"percent\":%d" */
void event(const char *name, const char *fields, ...) PRINTF_LIKE(2, 3);

/* report an error (human message, or an "error" event in JSON mode) and return code */
int fail(int code, const char *fmt, ...) PRINTF_LIKE(2, 3);

/* report a non-fatal warning (human message, or a "warning" event in JSON mode) */
void warn(const char *fmt, ...) PRINTF_LIKE(1, 2);

/* highlight escape codes, or "" when stderr is not a terminal that understands them */
const char *highlight(void);
const char *highlightEnd(void);

/* progress (fraction 0..1): a bar on stderr, or "progress" events in JSON mode;
   only reports when the whole-number percentage changes */
void progress(double fraction);

/* allocation wrappers: print an error and exit if memory runs out, so callers never see NULL */
void *xmalloc(size_t size);
void *xcalloc(size_t count, size_t size);
void *xrealloc(void *ptr, size_t size);

#endif
