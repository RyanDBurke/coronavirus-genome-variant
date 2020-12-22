all: run

run: build.c
	gcc -g -Wall -Wno-format -O2 -o run build.c ./utils/* -lz

clean:
	$(RM) run
