#-----------------------------------#
#	Makefile for GNU C++ Compiler	#
#-----------------------------------#
CCP = g++
CC = gcc
C_OMP = -fopenmp 
CFLAGS = -Wall -lgsl -lgslcblas -lm -O4 -fpic -funroll-loops -lgomp

FILES = metapop.c beetles.c crowley.c gillespie.c gillespie_recording.h gamma_beetles.c main.c
TRIBOL = tribolium.cpp kde.cpp 
WARNING = warning_signals.c odeintegrators.c gillespie.c correlation.c saddle_node.c record.c pars.c main.c 
CROWLEY = crowley.c gillespie.c main.c
BEETLES = beetles.c gillespie.c main.c
META = metapop.c gillespie.c gillespie_recording.c main.c  
SADDLENODE = saddle_node_direct.c gillespie.c gillespie_recording.c main.c



#build:
#	$(CC) $(CFLAGS) $(C_OMP) -o ibm.exe $(FILES)

saddlenode:
	$(CC) $(CFLAGS) $(C_OMP) -o saddlenode.exe $(SADDLENODE)

warning: 
	$(CC) $(CFLAGS) $(C_OMP) -o warning.exe $(WARNING)

crowley:
	$(CC) $(CFLAGS) $(C_OMP) -o crowley.exe $(CROWLEY)

beetles:
	$(CC) $(CFLAGS) $(C_OMP) -o beetles.exe $(BEETLES)

tribol:
	$(CCP) $(CFLAGS) $(C_OMP) -o tribol.exe $(TRIBOL)

#meta: 
#	$(CC) $(CFLAGS) $(C_OMP) -o meta.exe $(META)

clean:
	rm -f *.o 



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

