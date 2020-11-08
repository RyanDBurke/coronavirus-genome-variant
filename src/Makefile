all: fmmap

fmmap: fmmap.c
	gcc -g -Wall -O2 -o fmmap fmmap.c ./utils/align.c ./utils/fmIndex.c ./utils/misc.c -lz

clean:
	$(RM) fmmap
