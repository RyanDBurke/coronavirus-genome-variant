all: fmmap test z x

fmmap: fmmap.c
	gcc -g -Wall -o fmmap fmmap.c

test: test.c
	gcc -g -Wall -o test test.c

z: z.c
	gcc -g -Wall -o z z.c

x: x.c
	gcc -g -Wall -o x x.c

clean:
	$(RM) fmmap
	$(RM) test
	$(RM) z
	$(RM) x
