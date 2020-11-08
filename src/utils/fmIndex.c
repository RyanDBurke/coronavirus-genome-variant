#include "../lib/fmIndex.h"
KSEQ_INIT(gzFile, gzread)


int fmIndex(FM *fm, char *reference, char *output) {

    /* fasta stuff */
    gzFile fp;
	kseq_t *seq;
	int l;
	if (reference == NULL) {
		fprintf(stderr, "%s is NULL\n", reference);
		return 1;
	}
	fp = gzopen(reference, "r");
	seq = kseq_init(fp);

    /* parse .fa file contains our reference sequence and build fm-index */
    while ((l = kseq_read(seq)) >= 0) {

        /* name */
        fm->name = malloc(strlen(seq->name.s) + 1);
        strcpy(fm->name, "MN988713.1");

        /* seq */
        fm->seq = malloc(strlen(seq->seq.s) + 1);
        strcpy(fm->seq, seq->seq.s);
        
        /* length of seq */
        fm->length = strlen(seq->seq.s);

        /* suffix array */
        fm->suffixArray = (buildSuffixArray(seq->seq.s, strlen(seq->seq.s)));

        /* burrows-wheeler matrix (bwm) */
        fm->bwm = (buildBWM(seq->name.s, strlen(seq->seq.s)));

        /* burrows-wheeler transform (bwt) */
        fm->bwt = malloc(strlen(seq->seq.s) + 1);
        buildBWT(fm->bwt, fm->bwm, strlen(seq->seq.s));

        /* (F)irst and (L)ast columns of bwm */
        fm->F = malloc(strlen(seq->seq.s) + 1);
        fm->L = malloc(strlen(seq->seq.s) + 1);
        getFL(fm->F, fm->L, fm->bwm, strlen(seq->seq.s));

        /* Fast-rank arrays */
        /* to be implemented */

    }

    /* write FM-Index to output */
    FILE *f;
    int seqLimit = 50; /* adjust seq-length for writing to output to your liking */
    if (fm->length <= seqLimit) { 
        f = fopen(output, "w");
        if (f == NULL) {
            printf("error opening file: ");
            printf("\033[1;31m");
            printf("%s\n", output);
            printf("\033[0m");
            exit(1);
        }
        writeFM(fm, f, fm->length);
        printf("> You can find the serialized FM-Index for \"%s\" in: ", reference);
        printf("\033[1;31m");
        printf("%s\n\n", output);
        printf("\033[0m");
    } else {
        printf("\t> Sequences over %d in length are not written to file\n", seqLimit);
        printf("\t> Your sequence length: %d\n", fm->length);
        printf("\t> Use ");
        printf("\033[1;31m");
        printf("./fmmap default ");
        printf("\033[0m");
        printf("to see a serialized FM-Index.\n\n");
    }

    kseq_destroy(seq);
	gzclose(fp);
	return 0;
}

/*****************************/
/* AUX FUNTIONS FOR FM-INDEX */
/*****************************/

/* builds suffix array */
int *buildSuffixArray(char *seq, int length) {

    /* store suffixes and their offset */
    SA suffixes[length];

    /* starting adding SA objects to array */
    for (int i = 0; i < length; i++) {

        /* offset */
        suffixes[i].offset = i;

        /* build suffix */
        char* currentSuffix = (seq + i);
        suffixes[i].suffix = malloc(strlen(currentSuffix) + 1);
        strcpy(suffixes[i].suffix, currentSuffix);
    }

    /* sort suffixes */
    qsort(suffixes, length, sizeof(SA), &cmpSA);

    /* build suffix array */
    int *sa = malloc(sizeof(int) * length);
    for (int i = 0; i < length; i++) {

        /* store offset in our suffix array */
        sa[i] = suffixes[i].offset;
    }

    return sa;
}

/* build burrows-wheeler matrix */
char **buildBWM(char *seq, int length) {

    /* store rotations and their offset */
    R bwMatrix[length];

    /* starting adding rotations */
    for (int i = 0; i < length; i++) {

        /* offset */
        bwMatrix[i].offset = i;

        /* build rotation */
        char *prefix = malloc(i + 1);
        char *suffix = (seq + i);
        bwMatrix[i].rotation = malloc(length + 1);
        strcpy(bwMatrix[i].rotation, suffix);
        strncpy(prefix, seq, i);
        prefix[i] = '\0'; // add null terminator to prefix (crucial!)
        strcat(bwMatrix[i].rotation, prefix);

        /* free prefix, we don't need it anymore */
        free(prefix);
    }

    /* sort rotations */
    qsort(bwMatrix, length, sizeof(R), &cmpBMW);

    /* stores the actual matrix (rotations sorted) */
    char **matrix = malloc(sizeof(char *) * (length + 1));

    /* fill matrix */
    for (int i = 0; i < length; i++) {
        matrix[i] = malloc(sizeof(char) * length);
        strcpy(matrix[i], bwMatrix[i].rotation);
    }
    
    return matrix;
}

/* builds burrows-wheeler transform (i.e last column) */
void buildBWT(char *bwt, char **BWM, int length) {
    for (int i = 0; i < length; i++) {
        bwt[i] = BWM[i][length - 1];
    }
}

/* gets (F)irst and (L)ast columns of our burrows-wheeler matrix */
void getFL(char *F, char *L, char **BWM, int length) {
    for (int i = 0; i < length; i++) {
        F[i] = BWM[i][0];
        L[i] = BWM[i][length - 1];
    }
}

/* comparison sort function for suffix arrays */
int cmpSA(const void *a, const void *b) {
    const SA *da = (const SA *) a;
    const SA *db = (const SA *) b;

    return strcmp(da->suffix, db->suffix);
}

/* comparison sort function for bwm */
int cmpBMW(const void *a, const void *b) {
    const R *da = (const R *) a;
    const R *db = (const R *) b;

    return strcmp(da->rotation, db->rotation);
}

/* write our fm-index to output */
void writeFM(FM *fm, FILE *f, int length) {

    fprintf(f, "************\n");
    fprintf(f, "* FM-INDEX *\n");
    fprintf(f, "************\n\n");

    /* name */
    fprintf(f, "> %s\n", fm->name);

    /* sequence length */
    fprintf(f, "length: %d\n\n", fm->length);

    /* write sequence if its under 512 chars */
    fprintf(f, "Sequence\n");
    fprintf(f, "%s\n\n", fm->seq);

    /* suffix array */
    int lineCount = 0;
    fprintf(f, "Suffix Array\n");
    for (int i = 0; i < length; i++) {
        lineCount++;
        if (i == length - 1) {
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
    for (int i = 0; i < length; i++) {
        int offset = fm->suffixArray[i];
        char *suffix = (fm->seq + offset);

        if (i == length - 1) {
            fprintf(f, "%d:\t %s\n\n", i, suffix);
        } else {
            fprintf(f, "%d:\t %s\n", i, suffix);
        }
    }

    /* bwm */
    fprintf(f, "Burrows-Wheeler Matrix\n");
    for (int i = 0; i < length; i++) {
        if (i == length - 1) {
            fprintf(f, "%s\n\n", fm->bwm[i]);
        } else {
            fprintf(f, "%s\n", fm->bwm[i]);
        }
    }

    /* bwt */
    fprintf(f, "BWT: ");
    fprintf(f, "%s\n\n", fm->bwt);

    /* F and L */
    fprintf(f, "F: ");
    fprintf(f, "%s\n", fm->F);
    fprintf(f, "L: ");
    fprintf(f, "%s\n", fm->L);

    fclose(f);
}

/* prints suffix array */
void printSA(int *sa, int length) {

    printf("Suffix Array\n");
    for (int i = 0; i < length; i++) {
        if (i == length - 1) {
            printf("%d\n\n", sa[i]);
        } else {
            printf("%d, ", sa[i]);
        }
    }
}

/* prints bwm */
void printBWM(char **b, int length) {

    printf("Burrows-Wheeler Matrix\n");
    for (int i = 0; i < length; i++) {
        if (i == length - 1) {
            printf("%s\n\n", b[i]);
        } else {
            printf("%s\n", b[i]);
        }
    }
}

/* prints bwt */
void printBWT(char *bwt, int length) {
    printf("BWT: ");
    printf("%s\n\n", bwt);
}

/* prints F and L */
void printFL(char *F, char *L, int length) {
    printf("F: ");
    printf("%s\n", F);

    printf("L: ");
    printf("%s\n\n", L);
}

/* destroys our FM-Index*/
void destroy(FM *fm) {
    free(fm->bwt);
    free(fm->F);
    free(fm->L);
    free(fm->bwm);
    free(fm->suffixArray);
    free(fm->seq);
    free(fm->name);
    free(fm);
}
