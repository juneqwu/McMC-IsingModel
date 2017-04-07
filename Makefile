CC = gcc
CFLAGS = -fopenmp 
LDFLAGS = -lm
EXECS = ising

all: ${EXECS}

ising: ising.c
	${CC} ${CFLAGS} ising.c -o ising ${LDFLAGS} 

clean:

	rm -f ${EXECS} output00.txt output01.txt output02.txt output03.txt output04.txt output05.txt output06.txt output07.txt
