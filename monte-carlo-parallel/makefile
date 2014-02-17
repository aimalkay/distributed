CC = gcc
MPICC = mpicc
CFLAGS = -Wall -O2 -g -std=c99


OBJECTS2 = fcn.o monte-carlo.o main.o parallel-monte-carlo.o check-result.o

all: a1

a1: $(OBJECTS2)
	$(MPICC) -o a1 $(OBJECTS2) -lm 

monte-carlo.o: monte-carlo.c
	$(MPICC) $(CFLAGS) -c monte-carlo.c

fcn.o: fcn.c
	$(MPICC) $(CFLAGS) -c fcn.c

parallel-monte-carlo.o: parallel-monte-carlo.c
	$(MPICC) $(CFLAGS) -c parallel-monte-carlo.c

check-result.o: check-result.c
	$(MPICC) $(CFLAGS) -c check-result.c

main.o: main.c
	$(MPICC) $(CFLAGS) -c main.c

clean:
	rm $(OBJECTS1) $(OBJECTS2) a1