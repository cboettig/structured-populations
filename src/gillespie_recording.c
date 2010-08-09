/** record class for two-dimensional Gillespie simulations */
#include "gillespie.h"
#include "gillespie_recording.h"

record* record_alloc(size_t N, size_t replicates, double maxtime)
{
	record* myrecord = (record*) malloc(sizeof(record));
	myrecord->t_step = maxtime/N;
	myrecord->N = N;
	myrecord->s1 = (double*) calloc(replicates*(N+1), sizeof(double));
	myrecord->s2 = (double*) calloc(replicates*(N+1), sizeof(double));
	return myrecord;
}

void record_free(record * myrecord)
{
	free(myrecord->s1);
	free(myrecord->s2);
	free(myrecord);
}

/* A standard outcome function that just takes a flag to illustrate which indices of the vector (which pops) get birth or death events */ 
/* Doesn't follow the format of an event_fn object though, so would require an extra function call  */
double birth_out(void *ss, int i)
{ 
	double * s = (double *) ss;
	s[i] += 1;
	return 0;
}
double death_out(void * ss, int i)
{ 
	double * s = (double *) ss;
	s[i] -= 1;

	if(s[i] <= 0){
		fprintf(stderr, "extinction\n"); 
		return 1 ;
	}
	return 0 ; 
}


