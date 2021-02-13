all: run

run: ./src/main.c
	gcc -g -Wall -Wno-format -O2 -o run ./src/main.c ./utils/* -lz

clean:
	$(RM) run
