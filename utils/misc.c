#include "../lib/misc.h"

/* lowercase of a string */
char* lower(char* s) {
    for (int i = 0; i < strlen(s); i++) {
        s[i] = tolower(s[i]);
    }
    return s;
}

/* uppercase of a string */
char* upper(char* s) {
    for (int i = 0; i < strlen(s); i++) {
        s[i] = toupper(s[i]);
    }
    return s;
}

/* progress bar */
void pb(bool makeProgressBar, double percentageDone) {
    if(percentageDone >= 0 && percentageDone < .10) {
            printf("\rProgress: [#           ] %.0lf%%", percentageDone*100);
            fflush(stdout);
        } else if (percentageDone > .10 && percentageDone < .20) {
            printf("\rProgress: [##          ] %.0lf%%", percentageDone*100);
            fflush(stdout);
        } else if (percentageDone > .20 && percentageDone < .30) {
            printf("\rProgress: [###         ] %.0lf%%", percentageDone*100);
            fflush(stdout);
        } else if (percentageDone > .30 && percentageDone < .40) {
            printf("\rProgress: [####        ] %.0lf%%", percentageDone*100);
            fflush(stdout);
        } else if (percentageDone > .40 && percentageDone < .50) {
            printf("\rProgress: [#####       ] %.0lf%%", percentageDone*100);
            fflush(stdout);
        } else if (percentageDone > .50 && percentageDone < .60) {
            printf("\rProgress: [######      ] %.0lf%%", percentageDone*100);
            fflush(stdout);
        } else if (percentageDone > .60 && percentageDone < .70) {
            printf("\rProgress: [########    ] %.0lf%%", percentageDone*100);
            fflush(stdout);
        } else if (percentageDone > .70 && percentageDone < .80) {
            printf("\rProgress: [#########   ] %.0lf%%", percentageDone*100);
            fflush(stdout);
        } else if (percentageDone > .80 && percentageDone < .90) {
            printf("\rProgress: [##########  ] %.0lf%%", percentageDone*100);
            fflush(stdout);
        } else if (percentageDone > .90 && percentageDone < .98) {
            printf("\rProgress: [########### ] %.0lf%%", percentageDone*100);
            fflush(stdout);
        } else {
            printf("\rProgress: [############] %.0lf%%", percentageDone*100);
            fflush(stdout);
        }
}

/* print commands */
void commands() {

    /* Valid commands */
    printf("> Valid Commands:\n");
    printf("\033[1;31m");
    printf("   ./run\n");
    printf("   ./run covid 1K\n");
    printf("   ./run covid 10K\n");
    printf("   ./run covid 1M\n");
    printf("   ./run <reference.fa> <reads.fa[.gz]> <output.sam>\n\n");
    printf("\033[0m");
}

/* sets the paths for one of the "covid" presets */
static void covidPreset(char **ref_seq, char **reads, char *readsFile, double *READ, double count, char *label) {
    *ref_seq = "2019-nCoV.fa";
    *reads = readsFile;
    *READ = count;
    printf("[Aligning over %s reads]\n\n", label);
}

/* handle command-line args; returns 0 on success, non-zero for invalid usage */
int handleArgs(int argc, char **argv, char **ref_seq, char **indexOut, char **reads, char **alignOut, double *READ, bool *makeProgressBar) {

    /* defaults shared by every mode */
    *indexOut = "./fm-output/FMindex.txt";
    *alignOut = "./mappings/mapping.sam";
    *makeProgressBar = true;

    /* [./run] executes with small default inputs */
    if (argc == 1) {
        *ref_seq = "./default/ref-small.fa";
        *reads = "./default/reads-small.fa";
        *READ = 1;
        return 0;
    }

    /* [./run covid 1K|10K|1M] executes for coronavirus genome */
    if (argc == 3 && strcasecmp(argv[1], "covid") == 0) {
        if (strcasecmp(argv[2], "1K") == 0) {
            covidPreset(ref_seq, reads, "./reads/reads_1K.fa.gz", READ, 1000, "1,000");
            return 0;
        }
        if (strcasecmp(argv[2], "10K") == 0) {
            covidPreset(ref_seq, reads, "./reads/reads_10K.fa.gz", READ, 10000, "10,000");
            return 0;
        }
        if (strcasecmp(argv[2], "1M") == 0) {
            covidPreset(ref_seq, reads, "./reads/reads_1M.fa.gz", READ, 1000000, "1M");
            return 0;
        }
        return 1;
    }

    /* [./run <reference.fa> <reads.fa[.gz]> <output.sam>] explicit paths */
    if (argc == 4) {
        *ref_seq = argv[1];
        *reads = argv[2];
        *alignOut = argv[3];
        *READ = 0;                  /* unknown read count */
        *makeProgressBar = false;   /* so no progress bar */
        return 0;
    }

    return 1;
}
