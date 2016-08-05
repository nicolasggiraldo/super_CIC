CC = mpicc
CFLAGS = -Wall -I/home/$(USER)/local/include/ -O3
LFLAGS = -L/home/$(USER)/local/lib -lgsl -lgslcblas -lm -lfftw3
PROGRAM = super_CIC

$(PROGRAM):
	$(CC) -c $@.c $(CFLAGS)
	$(CC) $@.o $(LFLAGS) -o $@.x
	rm $@.o

clean:
	rm -rf *.out
	rm -rf *-
	rm -rf *.out
	rm -rf *#
	rm -rf *.o
	rm -rf *.a
	rm -rf *.so
	rm *.x
