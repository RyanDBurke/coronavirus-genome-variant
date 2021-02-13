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
    printf("   ./run covid 1M\n\n");
    printf("\033[0m");
    //printf("---------------------------------------------------------------------------------------------------\n\n");
}

/* handle command-line args */
void handleArgs(char **argv, char **ref_seq, char **indexOut, char **reads, char **alignOut, double *READ, bool *makeProgressBar) {

    /* [./run] executes with small default inputs */
    if (argv[1] == NULL) {
        *ref_seq = "./default/ref-small.fa";
        *indexOut = "./fm-output/FMindex.txt";
        *reads = "./default/reads-small.fa";
        *alignOut = "./mappings/mapping.sam";

        /* pass # of reads and make progress bar */
        *READ = 1;
        *makeProgressBar = true;
    }

    /* [./run covid 1k] executes for coronavirus genome with 1,000 reads*/
    else if (strcmp(lower(argv[1]), "covid") == 0 && (strcmp(argv[2], "1K") == 0 || strcmp(argv[2], "1k") == 0)) {
        *ref_seq = "2019-nCoV.fa";
        *indexOut = "./fm-output/FMindex.txt";
        *reads = "./reads/reads_1K.fa.gz";
        *alignOut = "./mappings/mapping.sam";

        printf("[Aligning over 1,000 reads]\n\n");

        /* pass # of reads and make progress bar */
        *READ = 1000;
        *makeProgressBar = true;
    }

    /* [./run covid 10k] executes for coronavirus genome with 10,000 reads*/
    else if (strcmp(lower(argv[1]), "covid") == 0 && (strcmp(argv[2], "10K") == 0 || strcmp(argv[2], "10k") == 0)) {
        *ref_seq = "2019-nCoV.fa";
        *indexOut = "./fm-output/FMindex.txt";
        *reads = "./reads/reads_10K.fa.gz";
        *alignOut = "./mappings/mapping.sam";

        printf("[Aligning over 10,000 reads]\n\n");

        /* pass # of reads and make progress bar */
        *READ = 10000;
        *makeProgressBar = true;
    }

    /* [./run covid 1M] executes for coronavirus genome with 1 Million reads */
    else if (strcmp(lower(argv[1]), "covid") == 0 && (strcmp(argv[2], "1M") == 0 || strcmp(argv[2], "1m") == 0)) {
        *ref_seq = "2019-nCoV.fa";
        *indexOut = "./fm-output/FMindex.txt";
        *reads = "./reads/reads_1M.fa.gz";
        *alignOut = "./mappings/mapping.sam";

        printf("[Aligning over 1M reads]\n\n");

        /* pass # of reads and make progress bar */
        *READ = 1000000;
        *makeProgressBar = true;
    } 
    
    /* just execute default*/
    else {
        *ref_seq = "./default/ref-small.fa";
        *indexOut = "./fm-output/FMindex.txt";
        *reads = "./default/reads-small.fa";
        *alignOut = "./mappings/mapping.sam";

        /* pass # of reads and make progress bar */
        *READ = 1;
        *makeProgressBar = true;
    }
}