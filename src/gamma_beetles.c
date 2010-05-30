#include "gillespie.h"


// K substages
// K egg transitions, K larva transitions, K pupa transitions,
// K egg deaths, K larva deaths, K pupa deaths, 1 adult death  starts at 3K+1 begins deaths d_start
// 1 birth -- stored in 6K+2
// Totals 6K+2 events
// 
// //
// X[3K+1] = K egg classes, K larva classes, K pupa classes, adult class, 
// Totals: 3K+1 states


//parameters = {ae, al, ap, ue, ul, up, ua, b, K, cle, cae, cap}

/** Calculate rates */
void rates_calc(double * rates, const int * states, const double * parameters){
	double transition_rate[3] = {parameters[0], parameters[1], parameters[2]};
	double death_rate[4] = {parameters[3], parameters[4], parameters[5], parameters[6]};
	double b = parameters[7];
	int K = parameters[8];
	const size_t n_rates = 6*K+2;
	const size_t n_states = 3*K+1;
	int i;

	/* sum classes for cannibalism */
	int classes[4] = {0,0,0,0};
	for(i = 0; i < n_states; i++)
	{
		classes[i / K] += states[i] ;
	}

	/* rates calc for eggs, larva, pupa */
	for(i = 0; i < n_states-1; i++){
		rates[i] = transition_rate[i/K] * states[i];	/* transitions: 0 : 3K-1 */
		rates[i+n_states-1] = death_rate[i/K] * states[i];/* mortality: 3K :6K-1 */
		if(i < K) 
		{	/* additional egg morality from cannibalism */
			rates[i+n_states-1] += parameters[9]*classes[1]*states[i]/parameters[12] +
				parameters[10]*classes[3]*states[i]/parameters[12];
		}
		if(i >= 2*K && i < 3*K)
		{	/* additional pupa morality from cannibalism */
			rates[i+n_states-1] += parameters[11]*classes[3]*states[i]/parameters[12];
		}
	}
	/* adults births */
	rates[n_rates-2] = death_rate[3]*states[n_states-1];
	rates[n_rates-1] = b*states[n_states-1];

	/* cumsum the rates */
	for(i = 1; i < n_rates; i++)
	{
		rates[i] += rates[i-1]; 
	}
}

/* event handling */
int outcome(int * states, const double * parameters, int event){
	int K = parameters[8];
	if (event < 3*K){ /* transitions */
		--states[event];
		++states[event+1];
	} else if (event < 6*K+1){ /* deaths */
		--states[event - 3*K ];
		if( states[event - 3*K] < 0) printf("whoops! %d %d\n", event, event - 3*K);
	} else if (event == 6*K+1){ /* birth */
		++states[0];
	} else { printf("error: event is: %d\n", event); }
	return 0; // switch to 1 to stop simulation
}

void 
fixed_interval_fn(const double t, const int * states, const double * parameters, void * my_record)
{
	double * rec = (double *) my_record;
	if( t > rec[0]){
		int classes[4] = {0,0,0,0};
		int K = parameters[8];
		const size_t n_states = 3*K+1;
		int i;
		for(i = 0; i < n_states; i++)
		{
			classes[i / K] += states[i] ;
		}
		printf("%g %d %d %d %d\n", t, classes[0], classes[1], classes[2], classes[3]);
		rec[0] += rec[1];
	}
};


/** */
void reset_state(int * states, const int * inits, const int n_states)
{
	int i;
	for(i = 0; i < n_states; i++)
	{
		states[i] = inits[i];
	}
}


void 
gillespie_sim(
	const int * inits, 
	const double * parameters, 
	const size_t n_rates, 
	const size_t n_states, 
	void * my_record, 
	const double max_time, 
	const size_t ensembles )
{
	/* intialize random number generator and seed with timestamp */
	gsl_rng *rng = gsl_rng_alloc (gsl_rng_default);
	#ifdef _RNDSEED
	gsl_rng_set(rng, time(NULL));
	#endif

	/* Gillespie simulation variables */
	int l,i,check;
	double lambda, t, tmp;

    /* Dynamically allocated private arrays must be declared inside 
	 * the parallel region.  */
	#pragma omp parallel shared(rng, inits, parameters, max_time, ensembles, my_record, n_rates, n_states) \
		private(lambda, t, tmp, i, check, l)
	{
		/* The vector to store cumulative sum of rates */
		double * rates_data = (double *) calloc (n_rates,sizeof(double) );
		int * states = (int *) calloc (n_states, sizeof(int) );
		/* Loop over ensembles, will be parallelized if compiled with -fopenmp */
		#pragma omp for
		for(l = 0; l < ensembles; l++){
			reset_state(states, inits, n_states);	
			t = 0; check=0;
			while(t < max_time){
				/* calculate time until next event */
			rates_calc(rates_data, states, parameters);
			lambda = rates_data[n_rates-1]; 
				t += gsl_ran_exponential(rng, 1/lambda);
			
				/* Events such as sampling occur at regular intervals */
				fixed_interval_fn(t, states, parameters, my_record);

				/* Execute the appropriate event */
				tmp = gsl_rng_uniform(rng);
				for(i=0;i<n_rates;i++){
					if( tmp < rates_data[i]/lambda ){
						check = outcome(states, parameters, i);
						break;
					} 
				}
			if(check) break;

		   }/* end evolution */
		}/* end ensembles */
		free(rates_data);
		free(states);
	} /* end parallel */
	gsl_rng_free(rng);
}



int main(void)
{
// Totals 6K+2 events
// X[3K+1] = K egg classes, K larva classes, K pupa classes, adult class, 
// Totals: 3K+1 states




	const int K = 10;
	const size_t n_rates = 6*K+2;
	const size_t n_states = 3*K+1;
	int * inits = (int *) calloc(n_states, sizeof(int));

	inits[0] = 100;
	/*                      {ae,  al,  ap, ue ,   ul,   up , ua , b, K,  cle,  cae,  cap, Vol*/
	double parameters[13] = {1.3, 0.1, 1.5, .00, .001, .00, .003, 5, K,  .2,   0.5, .100, 100};

	double max_time = 200;
	int ensembles = 1;

	double my_record[2] = {0, 10};
	
	double rates_data[n_rates];
	rates_calc(rates_data, inits, parameters);

//	fixed_interval_fn(0, inits, parameters, my_record);
//	outcome(inits, parameters, 6*K+1);

	gillespie_sim(inits, parameters, n_rates, n_states, my_record, max_time, ensembles);
	free(inits);
	return 0;
}



