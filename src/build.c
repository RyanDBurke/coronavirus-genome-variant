/***************/
/*    BUILD    */
/***************/

#include "./lib/align.h"
#include "./lib/fmIndex.h"
#include "./lib/misc.h"

/* keep a global note on # of reads */
double READ = 0;
bool makeProgressBar = false;

/* runs align and fmIndex */
int main(int argc, char **argv) {

    /* Valid commands */
    commands();

    /* string path to each relevant file */
    char *ref_seq = argv[1];
    char *indexOut = argv[2];
    char *reads = argv[3];
    char *alignOut = argv[4];

    /* [./fmmap default] executes with small default inputs */
    if (argv[1] == NULL) {
        ref_seq = "./default/ref-small.fa";
        indexOut = "./fm-output/FMindex.txt";
        reads = "./default/reads-small.fa";
        alignOut = "./mappings/mapping.sam";

        /* pass # of reads and make progress bar */
        READ = 1;
        makeProgressBar = true;
    }

    /* [./fmmap covid] executes for coronavirus genome with 1,000 reads*/
    else if (strcmp(lower(argv[1]), "covid") == 0 && (strcmp(argv[2], "1K") == 0 || strcmp(argv[2], "1k") == 0)) {
        ref_seq = "2019-nCoV.fa";
        indexOut = "./fm-output/FMindex.txt";
        reads = "./reads/reads_1K.fa.gz";
        alignOut = "./mappings/mapping.sam";

        printf("[Aligning over 1,000 reads]\n");
        printf("This takes roughly 30 seconds to execute\n\n");

        /* pass # of reads and make progress bar */
        READ = 1000;
        makeProgressBar = true;
    }

    /* [./fmmap covid] executes for coronavirus genome with 10,000 reads*/
    else if (strcmp(lower(argv[1]), "covid") == 0 && (strcmp(argv[2], "10K") == 0 || strcmp(argv[2], "10k") == 0)) {
        ref_seq = "2019-nCoV.fa";
        indexOut = "./fm-output/FMindex.txt";
        reads = "./reads/reads_10K.fa.gz";
        alignOut = "./mappings/mapping.sam";

        printf("[Aligning over 10,000 reads]\n");
        printf("This takes roughly 2 minutes to execute\n\n");

        /* pass # of reads and make progress bar */
        READ = 10000;
        makeProgressBar = true;
    }

    /* [./fmmap covid 1M] executes for coronavirus genome with 1 Million reads */
    else if (strcmp(lower(argv[1]), "covid") == 0 && (strcmp(argv[2], "1M") == 0 || strcmp(argv[2], "1m") == 0)) {
        ref_seq = "2019-nCoV.fa";
        indexOut = "./fm-output/FMindex.txt";
        reads = "./reads/reads_1M.fa.gz";
        alignOut = "./mappings/mapping.sam";

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
    align(fm, reads, alignOut, READ, makeProgressBar);

    /* free our fm-index */
    destroy(fm);

   return 0;
}
