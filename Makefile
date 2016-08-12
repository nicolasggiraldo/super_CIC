CC = mpicc
CFLAGS = -Wall -I/home/$(USER)/local/include/ -O3
LFLAGS = -L/home/$(USER)/local/lib -lgsl -lgslcblas -lm -lfftw3
PROGRAM = super_CIC

$(PROGRAM):
	$(CC) -c $@.c $(CFLAGS)
	$(CC) $@.o $(LFLAGS) -o $@.x
	rm $@.o

VEL:
	$(CC) -c $(PROGRAM).c -DVEL $(CFLAGS)
	$(CC) $(PROGRAM).o $(LFLAGS) -o $(PROGRAM).x
	rm $(PROGRAM).o

clean:
	rm -rf *.out
	rm -rf *-
	rm -rf *.out
	rm -rf *#
	rm -rf *.o
	rm -rf *.a
	rm -rf *.so
	rm *.x
