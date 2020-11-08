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
    printf("\t> ./fmmap covid 1K\n");
    printf("\t> ./fmmap covid 10K\n");
    printf("\t> ./fmmap covid 1M\n");
    printf("\t> ./fmmap default\n");
    printf("\t> ./fmmap <reference-sequence>.fa <output file>.txt <reads>.fa.gz <output file>.sam\n");
    printf("\033[0m");
    printf("---------------------------------------------------------------------------------------------------\n\n");
}