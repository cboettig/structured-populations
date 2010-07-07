#-----------------------------------#
#	Makefile for GNU C++ Compiler	#
#-----------------------------------#
CCP = g++
CC = gcc
C_OMP = -fopenmp 
CFLAGS = -Wall -lgsl -lgslcblas -lm -O3 -fpic -g
TRIBOL = tribolium.cpp kde.cpp 
WARNING = warning_signals.c odeintegrators.c gillespie.c correlation.c saddle_node.c record.c pars.c main.c 
CROWLEY = crowley.c gillespie.c main.c
BEETLES = beetles.c gillespie.c main.c
FILES = metapop.c beetles.c crowley.c gillespie.c hastings.c gamma_beetles.c main.c

build:
	$(CC) $(CFLAGS) $(C_OMP) -o ibm.exe $(FILES)

warning: 
	$(CC) $(CFLAGS) $(C_OMP) -o warning.exe $(WARNING)

crowley:
	$(CC) $(CFLAGS) $(C_OMP) -o crowley.exe $(CROWLEY)

beetles:
	$(CC) $(CFLAGS) $(C_OMP) -o beetles.exe $(BEETLES)

tribol:
	$(CCP) $(CFLAGS) $(C_OMP) -o tribol.exe $(TRIBOL)

clean:
	rm -f *.o *.so



### run make profile, then execute program, then run gprof -a $(EXE) gmon.out
#profile:  
#	$(CC) $(CFLAGS) -g -pg -o $(EXE) $(FILES)
#
#debug: 
#	$(CC) $(CFLAGS) -g -o $(EXE) $(FILES)
#
#parallel: $(WARNING)
#	$(CC) $(CFLAGS) $(C_OMP) -o $(EXE) $(WARNING) 
#
#

