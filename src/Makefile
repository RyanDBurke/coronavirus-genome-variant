all: fmmap test

fmmap: fmmap.c
	gcc -g -Wall -o fmmap fmmap.c

test: test.c
	gcc -g -Wall -o test test.c

clean:
	$(RM) fmmap
	$(RM) test
