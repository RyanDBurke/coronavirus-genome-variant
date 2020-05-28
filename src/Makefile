all: fmmap

fmmap: fmmap.c
	gcc -g -Wall -O2 -o fmmap fmmap.c -lz

clean:
	$(RM) fmmap
	rm -f *.o kseq_test
