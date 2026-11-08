CC      ?= gcc
CFLAGS  ?= -g -Wall -Wno-format -O2
LDLIBS   = -lz
SRCS     = utils/align.c utils/fmIndex.c utils/misc.c
HEADERS  = $(wildcard lib/*.h)
SANITIZE = -fsanitize=address,undefined -fno-omit-frame-pointer
PYTHON  ?= python3
OUT      = tests/out

all: run

run: src/main.c $(SRCS) $(HEADERS)
	$(CC) $(CFLAGS) -o $@ src/main.c $(SRCS) $(LDLIBS)

# unit tests run under AddressSanitizer + UndefinedBehaviorSanitizer
tests/test_units: tests/test_units.c $(SRCS) $(HEADERS)
	$(CC) -g -O1 -Wall -Wno-format $(SANITIZE) -o $@ tests/test_units.c $(SRCS) $(LDLIBS)

# regenerate the end-to-end fixture (output is committed; only needed if the generator changes)
fixture:
	$(PYTHON) tests/make_fixture.py 2019-nCoV.fa tests/fixtures/reads.fa tests/fixtures/truth.tsv

test: run tests/test_units
	ASAN_OPTIONS=detect_leaks=0 ./tests/test_units
	@mkdir -p $(OUT)
	./run 2019-nCoV.fa tests/fixtures/reads.fa $(OUT)/fixture.sam
	$(PYTHON) tests/validate_sam.py --ref 2019-nCoV.fa --sam $(OUT)/fixture.sam --truth tests/fixtures/truth.tsv
	./run 2019-nCoV.fa reads/reads_1K.fa.gz $(OUT)/reads_1K.sam
	$(PYTHON) tests/validate_sam.py --ref 2019-nCoV.fa --sam $(OUT)/reads_1K.sam --min-mapped $(MIN_MAPPED_1K)

# minimum fraction of the simulated 1K reads that must map
MIN_MAPPED_1K ?= 0.40

clean:
	$(RM) run tests/test_units
	$(RM) -r $(OUT)

.PHONY: all test fixture clean
