CC = mpicc
CFLAGS = -O2 -fopenmp
OBJECTS1 = main.o

all: multtable
multtable: $(OBJECTS1)
	$(CC) $(CFLAGS) -o $@ $? -lm
clean:
	@rm -rf $(OBJS) multtable *~ *.bak
