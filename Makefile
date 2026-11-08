CC       ?= gcc
CFLAGS   ?= -O2 -g
WARNINGS  = -Wall -Wextra -Werror
CPPFLAGS += -Iinclude -isystem include/vendor
DEPFLAGS  = -MMD -MP
LDLIBS    = -lz
SANITIZE  = -fsanitize=address,undefined -fno-omit-frame-pointer
PYTHON   ?= python3

BUILD     = build
OUT       = tests/out
LIB_SRCS  = src/align.c src/fmIndex.c src/misc.c
LIB_OBJS  = $(LIB_SRCS:src/%.c=$(BUILD)/%.o)
SAN_OBJS  = $(LIB_SRCS:src/%.c=$(BUILD)/asan/%.o)

# minimum fraction of the simulated 1K reads that must map
MIN_MAPPED_1K ?= 0.99

all: run

# ---- build ---------------------------------------------------------------

run: $(BUILD)/main.o $(LIB_OBJS)
	$(CC) $(CFLAGS) -o $@ $^ $(LDLIBS)

$(BUILD)/%.o: src/%.c | $(BUILD)
	$(CC) $(CPPFLAGS) $(DEPFLAGS) $(WARNINGS) $(CFLAGS) -c -o $@ $<

# AddressSanitizer + UndefinedBehaviorSanitizer builds, for the tests and `make asan`
$(BUILD)/asan/%.o: src/%.c | $(BUILD)/asan
	$(CC) $(CPPFLAGS) $(DEPFLAGS) $(WARNINGS) -g -O1 $(SANITIZE) -c -o $@ $<

$(BUILD)/asan/test_units.o: tests/test_units.c | $(BUILD)/asan
	$(CC) $(CPPFLAGS) $(DEPFLAGS) $(WARNINGS) -g -O1 $(SANITIZE) -c -o $@ $<

$(BUILD)/test_units: $(BUILD)/asan/test_units.o $(SAN_OBJS)
	$(CC) $(SANITIZE) -o $@ $^ $(LDLIBS)

$(BUILD)/run-asan: $(BUILD)/asan/main.o $(SAN_OBJS)
	$(CC) $(SANITIZE) -o $@ $^ $(LDLIBS)

asan: $(BUILD)/run-asan

$(BUILD) $(BUILD)/asan:
	mkdir -p $@

# ---- tests ---------------------------------------------------------------

test: run $(BUILD)/test_units $(BUILD)/run-asan
	./$(BUILD)/test_units
	@mkdir -p $(OUT)
	./$(BUILD)/run-asan --quiet -r 2019-nCoV.fa -q tests/fixtures/reads.fa -o $(OUT)/fixture.sam
	$(PYTHON) tests/validate_sam.py --ref 2019-nCoV.fa --sam $(OUT)/fixture.sam --truth tests/fixtures/truth.tsv
	./run --quiet -r 2019-nCoV.fa -q reads/reads_1K.fa.gz -o $(OUT)/reads_1K.sam
	$(PYTHON) tests/validate_sam.py --ref 2019-nCoV.fa --sam $(OUT)/reads_1K.sam --min-mapped $(MIN_MAPPED_1K)

# regenerate the end-to-end fixture (output is committed; only needed if the generator changes)
fixture:
	$(PYTHON) tests/make_fixture.py 2019-nCoV.fa tests/fixtures/reads.fa tests/fixtures/truth.tsv

clean:
	$(RM) run
	$(RM) -r $(BUILD) $(OUT)

.PHONY: all asan test fixture clean

-include $(wildcard $(BUILD)/*.d $(BUILD)/asan/*.d)
