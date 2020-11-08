all: run

run: fmmap.c
	gcc -g -Wall -O2 -o run fmmap.c ./utils/align.c ./utils/fmIndex.c ./utils/misc.c -lz

clean:
	$(RM) run
