#include "fmmap.h"
#include "kseq.h"
KSEQ_INIT(gzFile, gzread)

/* keep a global note on # of reads */
/* only works if a preset command is called */
double READ = 0;
bool makeProgressBar = false;

/* runs align and fmIndex */
int main(int argc, char **argv) {

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

    /* string path to each relevant file */
    char *ref_seq = argv[1];
    char *indexOut = argv[2];
    char *reads = argv[3];
    char *alignOut = argv[4];

    /* [./fmmap default] executes with small default inputs */
    if (strcmp(lower(argv[1]), "default") == 0) {
        ref_seq = "./Default/ref-small.fa";
        indexOut = "./FM-output/FMindex.txt";
        reads = "./Default/reads-small.fa";
        alignOut = "./Mappings/mapping.sam";

        /* pass # of reads and make progress bar */
        READ = 1;
        makeProgressBar = true;
    }

    /* [./fmmap covid] executes for coronavirus genome with 1,000 reads*/
    else if (strcmp(lower(argv[1]), "covid") == 0 && (strcmp(argv[2], "1K") == 0 || strcmp(argv[2], "1k") == 0)) {
        ref_seq = "2019-nCoV.fa";
        indexOut = "./FM-output/FMindex.txt";
        reads = "./Reads/reads_1K.fa.gz";
        alignOut = "./Mappings/mapping.sam";

        printf("[Aligning over 1,000 reads]\n");
        printf("This takes roughly 30 seconds to execute\n\n");

        /* pass # of reads and make progress bar */
        READ = 1000;
        makeProgressBar = true;
    }

    /* [./fmmap covid] executes for coronavirus genome with 10,000 reads*/
    else if (strcmp(lower(argv[1]), "covid") == 0 && (strcmp(argv[2], "10K") == 0 || strcmp(argv[2], "10k") == 0)) {
        ref_seq = "2019-nCoV.fa";
        indexOut = "./FM-output/FMindex.txt";
        reads = "./Reads/reads_10K.fa.gz";
        alignOut = "./Mappings/mapping.sam";

        printf("[Aligning over 10,000 reads]\n");
        printf("This takes roughly 2 minutes to execute\n\n");

        /* pass # of reads and make progress bar */
        READ = 10000;
        makeProgressBar = true;
    }

    /* [./fmmap covid 1M] executes for coronavirus genome with 1 Million reads */
    else if (strcmp(lower(argv[1]), "covid") == 0 && (strcmp(argv[2], "1M") == 0 || strcmp(argv[2], "1m") == 0)) {
        ref_seq = "2019-nCoV.fa";
        indexOut = "./FM-output/FMindex.txt";
        reads = "./Reads/reads_1M.fa.gz";
        alignOut = "./Mappings/mapping.sam";

        printf("[Aligning over 1M reads]\n");
        printf("This takes roughly ~2.5 Hours to execute\n\n");

        /* pass # of reads and make progress bar */
        READ = 1000000;
        makeProgressBar = true;
    }

    /* our FM-Index struct */
    FM *fm = malloc(sizeof(FM));

    /* run fmIndex */
    fmIndex(fm, ref_seq, indexOut);

    /* run aligner */
    align(fm, reads, alignOut);

    /* free our fm-index */
    destroy(fm);

   return 0;
}

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
        printf("to see what a serialized FM-Index looks like.\n");
        printf("\t> If you really want to see FM-Index for sequences over %d in length you can adjust it on line 151 in the file ", seqLimit);
        printf("\033[1;31m");
        printf("fmmap.c\n\n");
        printf("\033[0m");
    }

    kseq_destroy(seq);
	gzclose(fp);
	return 0;
}


int align(FM *fm, char *reads, char *output) {

    /* keep track of our progress */
    double progress = 0.0;

    /* fasta stuff */
    gzFile fp;
	kseq_t *seq;
	int l;
	if (reads == NULL) {
		fprintf(stderr, "%s is NULL\n", reads);
		return 1;
	}
	fp = gzopen(reads, "r");
	seq = kseq_init(fp);

    double ninf = -INFINITY;
    int gap = -5; /* gap penalty */

    /* open our output file (.sam) */
    FILE *f;
    f = fopen(output, "w");
    f = fopen(output, "a");

    /* .sam header */
    fprintf(f, "@HD	VN:1.0	SO:unsorted\n");
    fprintf(f, "@SQ	SN:%s	LN:%d\n", fm->name, fm->length);

    /* parse .fa file containing our n-amount of 100bp reads */
    while ((l = kseq_read(seq)) >= 0) {

        /* progress bar */
        progress++;
        double percentageDone = progress / READ;
        pb(true, percentageDone);

        /* an array where we'll store our alignments */
        Alignment alignments[strlen(seq->seq.s)];

        double bestScore = ninf;
        int skip = seedSkip(strlen(seq->seq.s));
        int alignmentLength = -1;

        /* for each (read-length / 5.0) seed (20bp seed in our case, since we have 100bp reads) */
        for (int seedStart = 0; seedStart < strlen(seq->seq.s); seedStart += skip) {
            int seedEnd = min(strlen(seq->seq.s), seedStart + skip);

            /* finding match interval */
            Interval *interval = malloc(sizeof(Interval));
            int matchLength = 0;
            char *seed = malloc(strlen(seq->seq.s) + 1);
            substring(seed, seq->seq.s, seedStart, seedEnd);
            getInterval(interval, &matchLength, fm, seed);

            /* if there wasn't a match, break out of this seed */
            if (matchLength == 0) {
                break;
            }

            /* reference position */
            int *refPos = malloc((interval->end - interval->start));
            int refPosLength = referencePos(refPos, interval, matchLength, fm, seedEnd);            

            /* fitting alignment: add it to alignments array and return that array length */
            alignmentLength = alignment(alignments, seq->seq.s, fm->seq, refPos, refPosLength, gap, bestScore);

            /* free memory */
            free(refPos);
            free(seed);
            free(interval);
        }

        /* for each alignment in alignments, write to .sam file */
        for (int a = 0; a < alignmentLength; a++) {
            if (f == NULL) {
                printf("error opening file: ");
                printf("\033[1;31m");
                printf("%s\n", output);
                printf("\033[0m");
                exit(1);
            } else {

                /* our SAM struct we will fill */
                SAM *s = malloc(sizeof(SAM));

                /* QNAME */
                s->QNAME = malloc(strlen(seq->name.s) + 1);
                strcpy(s->QNAME, seq->name.s);

                /* FLAG */
                s->FLAG = 0;

                /* RNAME */
                s->RNAME = malloc(strlen(fm->name) + 1);
                strcpy(s->RNAME, fm->name);

                /* POS */
                s->POS = alignments[a].pos;

                /* MAPQ */
                s->MAPQ = 255;

                /* CIGAR */
                s->CIGAR = malloc(strlen(alignments[a].cigar) + 1);
                strcpy(s->CIGAR, alignments[a].cigar);

                /* RNEXT */
                s->RNEXT = malloc(2);
                strcpy(s->RNEXT, "*");

                /* PNEXT */
                s->PNEXT = 0;

                /* TLEN */
                s->TLEN = strlen(seq->seq.s);

                /* SEQ */
                s->SEQ = malloc(strlen(seq->seq.s) + 1);
                strcpy(s->SEQ, seq->seq.s);

                /* QUAL */
                s->QUAL = malloc(2);
                strcpy(s->QUAL, "*");
                
                /* write SAM */
                writeSAM(s, f);
                destroySAM(s);
            }
        }      
    }

    fclose(f);

    printf("\n\n> You can find the SAM-formatted alignments in: ");
    printf("\033[1;31m");
    printf("%s\n", output);
    printf("\033[0m");

    /* successful */
    kseq_destroy(seq);
	gzclose(fp);
	return 0;
}

/**************************/
/* AUX FUNTIONS FOR ALIGN */
/**************************/

/* preforms backwards search and return an interval and matchLength */
void getInterval(Interval *interval, int *matchLength, FM *fm, char* seed) {

    /* reset interval to -1 */
    interval->start = -1;
    interval->end = -1;

    /* length of our original reference sequence (coronavirus genome) */
    int refSeqLength = fm->length;

    /* iterate thru suffix array to find matches to our seed */
    for (int i = 0; i < refSeqLength; i++) {

        /* index in coronavirus genome where this particular suffix begins */
        int offset = fm->suffixArray[i];

        /* current coronavirus suffix */
        char *suffix = (fm->seq + offset);

        /* substring of suffix */
        char *subSuffix = malloc(strlen(suffix) + 1);
        substring(subSuffix, suffix, 0, strlen(seed));

        /* match seed to suffix */
        int match = strcmp(subSuffix, seed);

        /* if they match and this is the first match */
        if (match == 0 && (interval->start == -1)) {
            interval->start = i;
        }

        /* if they don't match, but we've seen a match before -- then this is the end of interval */
        if (match != 0 && (interval->start != -1)) {
            interval->end = i;
            break; /* break loop, we have our interval! */
        }

        /* free subSuffix substring memory */
        free(subSuffix);
    }

    /* matchLength = 0 if no match found, otherwise matchLength = strlen(seed) */
    if ((interval->start == -1) && (interval->end == -1)) {
        *matchLength = 0;
    } else {
        *matchLength = strlen(seed);
    }
}

/* returns reference positions of current seed */
int referencePos(int *refPos, Interval *interval, int matchLength, FM *fm, int seedEnd) {

    /* interval: [start, end) */
    int start = interval->start;
    int end = interval->end;

    /* find reference positions in reference genome */
    int refPosLength = 0;
    for (int i = start; i < end; i++) {
        int pos = (fm->suffixArray[i]) - (seedEnd - matchLength);
        refPos[refPosLength] = pos;
        refPosLength++;
    }

    /* we're returning the length of refPos int-array */
    return refPosLength;
}

/* return a single alignment struct to output parameter A */
int alignment(Alignment alignments[], char *read, char *ref, int *refPos, int refPosLength, int gap, double bestScore) {

    /* X-axis is our coronavirus reference genome slice */
    /* Y-axis is our read we wish to align to our reference */

    /* index for our alignments-array */
    int alignmentIndex = 0;

    /* for each reference positions in our reference genome */
    for (int pos = 0; pos < refPosLength; pos++) {

        /* x and y-axis strings we will align */
        char *x = malloc(strlen(read) + 10 + 1); // len(read) + (2 * gap) + null terminator = 111
                                                 // (100)     + (10)      + (1)             = 111

        /* slice out reference genome with 5-gap on each side (110 total length) */
        /* we use max and min here to avoid going below 0 or above length in index */
        /* given that, our slice from reference genome is: 105 <= length <= 110 */
        substring(x, ref, max(refPos[pos] - 5, 0), min(refPos[pos] + strlen(read) + 5, strlen(ref)));

        /* setting y-axis is trivial */
        char *y = read;

        /* n and m are lengths of our current sliced-genome and read */
        int n = strlen(x); /* 105 <= length <= 110 */
        int m = strlen(read); /* length = 100 */

        /* our OPT-matrix containing edit-distance between x and y */
        int OPT[MAXROW][MAXCOL];

        /* our edit-distance between strings x and y */
        int score = editDistance(OPT, x, y, n, m, gap);

        /* only create an alignment for those with >= best score */
        if (score > bestScore) {

            /* update the best score */
            bestScore = score;

            /* we set this to zero to simulate "clearing" the array */
            alignmentIndex = 0;

            /* backtrace OPT-matrix to find alignment and return CIGAR string */
            int offset = -1;
            char *cigar = buildCigar(OPT, n, m, gap, x, y, &offset);

            /* add to alignments array */
            alignments[alignmentIndex].score = score;
            alignments[alignmentIndex].pos = (refPos[pos] + offset == 0) ? 1 : refPos[pos] + offset;
            // alignments[alignmentIndex].pos = refPos[pos] + offset;
            alignments[alignmentIndex].cigar = malloc(strlen(cigar) + 1);
            strcpy(alignments[alignmentIndex].cigar, cigar);
            alignmentIndex++;

        } else if (score == bestScore) { /* if its equal we add it to alignments-array */

            /* backtrace OPT-matrix to find alignment and return CIGAR string */
            int offset = -1;
            char *cigar = buildCigar(OPT, n, m, gap, x, y, &offset);

            /* add to alignments array */
            alignments[alignmentIndex].score = score;
            alignments[alignmentIndex].pos = (refPos[pos] + offset == 0) ? 1 : refPos[pos] + offset;
            // alignments[alignmentIndex].pos = refPos[pos] + offset;
            alignments[alignmentIndex].cigar = malloc(strlen(cigar) + 1);
            strcpy(alignments[alignmentIndex].cigar, cigar);
            alignmentIndex++;
        }

        /* free variables */
        free(x);
    }

    return alignmentIndex;
}

/* builds OPT-matrix, and returns OPT[n][m] -- which is the edit distance between x and y */
int editDistance(int OPT[MAXROW][MAXCOL], char *x, char *y, int n, int m, int gap) {

    /* (0, 0) is always 0 */
    OPT[0][0] = 0;

    /* add our initial gap penalties to the first column of each row */
    /* allow cost-free gap penalties */
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

            /* fil matrix entry */
            OPT[i][j] = maxAlign(
                    (score(x[i - 1], y[j - 1], gap) + OPT[i - 1][j - 1]),
                    (gap + OPT[i - 1][j]),
                    (gap + OPT[i][j - 1])
                );
        }
    }

    return OPT[n][m];
}

/* returns CIGAR-string */
char *buildCigar(int OPT[MAXROW][MAXCOL], int n, int m, int gap, char *x, char *y, int *offset) {

    /* find out where we start our backtrace, and set n */
    int temp = n;
    int maxScore = -1;
    for (int i = 0; i <= temp; i++) {
        if (OPT[i][m] >= maxScore) {
            maxScore = OPT[i][m];
            n = i;
        }
    }

    /* result array we'll return */
    int cigarIndex = 0;

    char *cigarTemp = malloc(n * m);
    while (true) {


        /* for fitting alignment it ends when we reach row = 0 */
        /* in our case, when m = 0 */
        if (m == 0) {

            /* Crucial
                * this decides where exactly our read aligns best to our reference genome
                * Allowing cost-free gaps at the beginning of our x-axis means our read
                    can align (technically) anywhere on on the slice of our reference genome
                    although it would typically be around of reference position
             */

            *offset = n;
            break;
        }

        /* find scores for neighboring entries */
        int up = gap + OPT[n - 1][m];
        int diagonal = score(x[n - 1], y[m - 1], gap) + OPT[n - 1][m - 1];
        int left = gap + OPT[n][m - 1];

        /* find max of those neighboring entries to find path*/
        int maxA = maxAlign(up, diagonal, left);

        /* build cigar string based on path () */
        if (maxA == diagonal) {
            n = n - 1;
            m = m - 1;
            cigarTemp[cigarIndex] = 'M';
            cigarIndex++;
        } else if (maxA == left) {
            m = m - 1;
            cigarTemp[cigarIndex] = 'I';
            cigarIndex++;
        } else if (maxA == up) {
            n = n - 1;
            cigarTemp[cigarIndex] = 'D';
            cigarIndex++;
        } else {
            printf("Error in buildCigar()\n");
            exit(1);
        }

    }

    /* null terminator */
    cigarTemp[cigarIndex] = '\0';

    /* now that we (almost) have our CIGAR-string, we need to format it properly */
    char *cigar = formatCigar(cigarTemp, strlen(cigarTemp) * 2 + 1);

    /* free used memory */
    free(cigarTemp);

    /* our complete cigar string */
    return cigar;
}

/* return score between two characters */
int score(char a, char b, int gap) {
    if (a == 'N' || b == 'N') { /* if you see 'N', assume it matches */
        return 0;
    } else if (a == b) {
        return 1;
    } else if (a != b) {
        return -1;
    } else if (a == '-' || b == '-') {
        return gap;
    } else {
        return 0;
    }
}

/* return substring of string from [start, end) */
void substring(char* result, char* string, int start, int end) {
    char *temp = string + start;
    int n = end - start;
    result = strncpy(result, temp, n);
    result[n] = '\0';
}

/* return correctly formatted CIGAR-string (i.e 3M2D2M) */
char *formatCigar(char *cigar, int length) {

    /* just make sure cigar input is uppercase */
    cigar = upper(cigar);

    /* our properly formatted CIGAR-string we will return */
    char *result = malloc(length + 1);
    int resultIndex = 0;

    /* build cigar string */
    int count = 0;
    char currentChar = '%'; /* random char, shouldn't matter */
    for (int i = 0; i < length; i++) {
        currentChar = cigar[i];

        for (int j = i; j < length; j++) {
            if (cigar[j] != currentChar) {
                if (count > 9) {
                    char array[64];
                    int myInteger = count;
                    sprintf(array, "%d", myInteger );

                    for (int k = 0; k < strlen(array); k++) {
                        result[resultIndex] = array[k];
                        resultIndex++;
                    }

                } else {
                    result[resultIndex] = count + '0';
                    resultIndex++;
                }
                
                result[resultIndex] = currentChar;
                resultIndex++;

                if (cigar[j + 1] == '\0') {
                    j = length;
                }
                break;
            } else {
                count++;

                if (cigar[j + 1] == '\0') {
                    if (count > 9) {
                        char array[64];
                        int myInteger = count;
                        sprintf( array, "%d", myInteger );

                        for (int k = 0; k < strlen(array); k++) {
                            result[resultIndex] = array[k];
                            resultIndex++;
                        }

                    } else {
                        result[resultIndex] = count + '0';
                        resultIndex++;
                    }

                    result[resultIndex] = currentChar;
                    resultIndex++;
                    j = length;
                }
            }
            i = j;
        }

        count = 0;
    }

    /* null-terminator */
    result[resultIndex] = '\0';    

    return result;
}

/* seed skip */
int seedSkip(int L) {
    int s = L / 5.0;
    return s;
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

/* print our OPT-matrix */
/* not useful for large genomes */
void printMatrix(int OPT[MAXROW][MAXCOL], int n, int m, char *x, char *y) {

    /* Y-axis */
    for (int i = 0; i <= m + 1; i++) {
        if (i == 0) {
            printf("\t");
        } else if (i == 1) {
            printf("Y-axis\t");
        } else {
            printf("%c\t", y[i - 2]);
        }
    }

    /* formatting */
    printf("\n");

    /* print our X-axis string and the OPT-matrix entries */
    for (int i = 0; i <= n; i++) {     
        if (i > 0) {
            printf("%c\t", x[i - 1]);
        } else {
            printf("X-axis\t");
        }
        for (int j = 0; j <= m; j++) {
            printf("%d\t", OPT[i][j]);
        }
        printf("\n");
    }
}

/* write to our .SAM file */
void writeSAM(SAM *sam, FILE *f) {
    fprintf(f, "%s\t", sam->QNAME);
    fprintf(f, "%d\t", sam->FLAG);
    fprintf(f, "%s\t", sam->RNAME);
    fprintf(f, "%d\t", sam->POS);
    fprintf(f, "%d\t", sam->MAPQ);
    fprintf(f, "%s\t", sam->CIGAR);
    fprintf(f, "%s\t", sam->RNEXT);
    fprintf(f, "%d\t", sam->PNEXT);
    fprintf(f, "%d\t", sam->TLEN);
    fprintf(f, "%s\t", sam->SEQ);
    fprintf(f, "%s\n", sam->QUAL);
}

/* free our SAM struct */
void destroySAM(SAM *sam) {

    /* some SAM fields are free'd elsewhere */
    // free(sam->QNAME);
    // free(sam->RNAME);
    free(sam->CIGAR);
    free(sam->RNEXT);
    // free(sam->PNEXT);
    // free(sam->SEQ);
    free(sam->QUAL);
    free(sam);
}

/* destroy our alignment array */
void destroyAlignment(Alignment *A) {}

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

/********/
/* MISC */
/********/

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
