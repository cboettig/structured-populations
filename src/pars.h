#include <stdlib.h>
#define N_EVENT_TYPES 2

/** Event functions must use this custom structure to represent 
 *  both the state of the system and any model parameters */
typedef struct {
	int n;
	int K;
	double e;
	double a;
	double h;
	int time_index;
	int hist_index;
	double checkpts[N_EVENT_TYPES];
	double start_polluting;
	double pollute_rate;
	double pollute_increment;
} pars;

void * pars_alloc(int N, int K, double e, double a, double h, double start_polluting, double pollute_rate, double pollute_increment);


