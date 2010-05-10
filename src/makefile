#-----------------------------------#
#	Makefile for GNU C++ Compiler	#
#-----------------------------------#
CC = gcc
C_OMP = -fopenmp 
CFLAGS = -Wall -lgsl -lgslcblas -lm -O3 -fpic
EXE = warning_signals.exe
FILES = warning_signals.c odeintegrators.c gillespie.c correlation.c main.c 

build: $(FILES)
	$(CC) $(CFLAGS) -o $(EXE) $(FILES)

parallel: $(FILES)
	$(CC) $(CFLAGS) $(C_OMP) -o $(EXE) $(FILES) 

profile: $(FILES)
	$(CC) $(CFLAGS) -g -pg -o $(EXE) $(FILES)

debug:  $(FILES)
	$(CC) $(CFLAGS) -g -o $(EXE) $(FILES)
	
