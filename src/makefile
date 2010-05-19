#-----------------------------------#
#	Makefile for GNU C++ Compiler	#
#-----------------------------------#
CCP = g++
CC = gcc
C_OMP = -fopenmp 
CFLAGS = -Wall -lgsl -lgslcblas -lm -O3 -fpic
BEETLES = tribolium.cpp kde.cpp 
WARNING = warning_signals.c odeintegrators.c gillespie.c correlation.c  record.c pars.c main.c 
CROWLEY = crowley.c gillespie.c 

build: 
	$(CC) $(CFLAGS) -o warning.exe $(WARNING)
	$(CC) $(CFLAGS) -o crowley.exe $(CROWLEY)
	$(CC) $(CFLAGS) -o beetles.exe $(BEETLES)

warning: $(WARNING)
	$(CC) $(CFLAGS) -o warning.exe $(WARNING)

crowley: $(CROWLEY)
	$(CC) $(CFLAGS) -o crowley.exe $(CROWLEY)

beetles: $(BEETLES)
	$(CC) $(CFLAGS) -o beetles.exe $(BEETLES)

clean:
	rm -f *.o *.so

# run make profile, then execute program, then run gprof -a $(EXE) gmon.out
profile:  
	$(CC) $(CFLAGS) -g -pg -o $(EXE) $(FILES)

debug: 
	$(CC) $(CFLAGS) -g -o $(EXE) $(FILES)

parallel: $(WARNING)
	$(CC) $(CFLAGS) $(C_OMP) -o $(EXE) $(WARNING) 



