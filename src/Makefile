all: fmmap test z

fmmap: fmmap.c
	gcc -g -Wall -o fmmap fmmap.c

test: test.c
	gcc -g -Wall -o test test.c

z: z.c
	gcc -g -Wall -o z z.c

clean:
	$(RM) fmmap
	$(RM) test
	$(RM) z
