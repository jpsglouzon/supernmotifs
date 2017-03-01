CC = g++

CFLAGS=-Wall -fexceptions -O3 -std=c++0x -m64

CFLAGS_o=-s -m64 -static -static-libgcc -static-libstdc++ 

supernmotifs: main.o RepMotif.o RepSupermotif.o RepWeightedMotif.o
	$(CC) -L. -o supernmotifs main.o RepMotif.o RepSupermotif.o RepWeightedMotif.o $(CFLAGS_o)

main.o: main.cpp
	$(CC) $(CFLAGS) -I. -c main.cpp -o main.o

RepMotif.o: RepMotif.cpp
	$(CC) $(CFLAGS) -I. -c RepMotif.cpp -o RepMotif.o

RepSupermotif.o: RepSupermotif.cpp
	$(CC) $(CFLAGS) -I. -c RepSupermotif.cpp -o RepSupermotif.o

RepWeightedMotif.o: RepWeightedMotif.cpp
	$(CC) $(CFLAGS) -I. -c RepWeightedMotif.cpp -o RepWeightedMotif.o

clean:
	rm -rf *.o


 



