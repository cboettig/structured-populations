#-----------------------------------#
#	Makefile for GNU C++ Compiler	#
#-----------------------------------#
CC = g++
CFLAGS = -Wall -lgsl -lgslcblas -lm
EXE = tribolium 
FILES = tribolium.cpp kde.c 

build: $(FILES)
	$(CC) $(CFLAGS) -o $(EXE) $(FILES)

clean:
	rm -f *.o *.so

# run make profile, then execute program, then run gprof -a $(EXE) gmon.out
profile:  $(FILES)
	$(CC) $(CFLAGS) -g -pg -o $(EXE) $(FILES)

debug: $(FILES)
	$(CC) $(CFLAGS) -g -o $(EXE) $(FILES)

