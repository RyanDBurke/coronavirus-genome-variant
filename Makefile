# GNU make's built-in default for CC is "cc"; use gcc unless CC is set explicitly
ifeq ($(origin CC),default)
CC        = gcc
endif
CFLAGS   ?= -O2 -g
WARNINGS  = -Wall -Wextra -Werror
CPPFLAGS += -Iinclude -isystem include/vendor
DEPFLAGS  = -MMD -MP
LDLIBS    = -lz
PYTHON   ?= python3

BUILD     = build
OUT       = tests/out

ifeq ($(OS),Windows_NT)
# native Windows (MSYS2 UCRT64): static run.exe with an embedded UTF-8 manifest; MinGW has no ASan/UBSan
EXE       = .exe
SANITIZE ?=
LDFLAGS  += -static
WINDRES  ?= windres
RES       = $(BUILD)/run.res
else
EXE       =
SANITIZE ?= -fsanitize=address,undefined -fno-omit-frame-pointer
RES       =
endif

LIB_SRCS  = src/align.c src/fmIndex.c src/misc.c
LIB_OBJS  = $(LIB_SRCS:src/%.c=$(BUILD)/%.o)
TEST_OBJS = $(LIB_SRCS:src/%.c=$(BUILD)/test/%.o)

# tests run the sanitizer build where there is one
TEST_RUN  = $(if $(strip $(SANITIZE)),$(BUILD)/run-asan$(EXE),./run$(EXE))

# minimum fraction of the simulated 1K reads that must map
MIN_MAPPED_1K ?= 0.99

all: run$(EXE)

# ---- build ---------------------------------------------------------------

run$(EXE): $(BUILD)/main.o $(LIB_OBJS) $(RES)
	$(CC) $(CFLAGS) $(LDFLAGS) -o $@ $^ $(LDLIBS)

$(BUILD)/%.o: src/%.c | $(BUILD)
	$(CC) $(CPPFLAGS) $(DEPFLAGS) $(WARNINGS) $(CFLAGS) -c -o $@ $<

$(BUILD)/run.res: win/run.rc win/run.manifest | $(BUILD)
	$(WINDRES) --include-dir win -O coff -o $@ $<

# test builds: AddressSanitizer + UndefinedBehaviorSanitizer where available
$(BUILD)/test/%.o: src/%.c | $(BUILD)/test
	$(CC) $(CPPFLAGS) $(DEPFLAGS) $(WARNINGS) -g -O1 $(SANITIZE) -c -o $@ $<

$(BUILD)/test/test_units.o: tests/test_units.c | $(BUILD)/test
	$(CC) $(CPPFLAGS) $(DEPFLAGS) $(WARNINGS) -g -O1 $(SANITIZE) -c -o $@ $<

$(BUILD)/test_units$(EXE): $(BUILD)/test/test_units.o $(TEST_OBJS)
	$(CC) $(SANITIZE) $(LDFLAGS) -o $@ $^ $(LDLIBS)

$(BUILD)/run-asan$(EXE): $(BUILD)/test/main.o $(TEST_OBJS)
	$(CC) $(SANITIZE) $(LDFLAGS) -o $@ $^ $(LDLIBS)

asan: $(BUILD)/run-asan$(EXE)

$(BUILD) $(BUILD)/test:
	mkdir -p $@

# ---- tests ---------------------------------------------------------------

test: run$(EXE) $(BUILD)/test_units$(EXE) $(TEST_RUN)
	./$(BUILD)/test_units$(EXE)
	@mkdir -p $(OUT)
	$(TEST_RUN) --quiet -r 2019-nCoV.fa -q tests/fixtures/reads.fa -o $(OUT)/fixture.sam
	$(PYTHON) tests/validate_sam.py --ref 2019-nCoV.fa --sam $(OUT)/fixture.sam --truth tests/fixtures/truth.tsv
	./run$(EXE) --quiet -r 2019-nCoV.fa -q reads/reads_1K.fa.gz -o $(OUT)/reads_1K.sam
	$(PYTHON) tests/validate_sam.py --ref 2019-nCoV.fa --sam $(OUT)/reads_1K.sam --min-mapped $(MIN_MAPPED_1K)
	$(PYTHON) tests/check_cli.py --run ./run$(EXE) --out $(OUT)

# regenerate the end-to-end fixture (output is committed; only needed if the generator changes)
fixture:
	$(PYTHON) tests/make_fixture.py 2019-nCoV.fa tests/fixtures/reads.fa tests/fixtures/truth.tsv

clean:
	$(RM) run run.exe
	$(RM) -r $(BUILD) $(OUT)

.PHONY: all asan test fixture clean

-include $(wildcard $(BUILD)/*.d $(BUILD)/test/*.d)
