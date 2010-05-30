#include "gillespie.h"


// K substages
// K egg transitions, K larva transitions, K pupa transitions,
// K egg deaths, K larva deaths, K pupa deaths, 1 adult death  starts at 3K+1 begins deaths d_start
// 1 birth -- stored in 6K+2
// rates[6K+2]
// 
// //
// X[3K+1] = K egg classes, K larva classes, K pupa classes, adult class




int K  = 10;
int N = 6*K+2;
d_start = 3*K+1;


parameters = {ae, al, ap, ue, ul, up, ua, b}

/** Calculate rates */
void rates_calc(double * rates, const double * states, const double * parameters){

	double transition_rate[3] = {parameters[0], parameters[1], parameters[2]};
	double death_rate[4] = {parameters[3], parameters[4], parameters[5], parameters[6]};
	double b = parameters[7]; 

	int i;
	for(i = 0; i < N-1; i++){
		rates[i] = transition_rate[i%K] * states[i];
		rates[i+d_start] = death_rate[i%K] * states[i];
	}
	rates[N] = b*states[N-1];

	/* cumsum the rates */
	for(i = 1; i < N; i++)
	{
		rates[i] += rates[i-1]; 
	}


}

/* event handling */
void outcome(double * state, int event){
	if (event < 3*K){
		--state[event];
		++state[event+1]; 
	} else if (event < 6*K){
		--state[event];
	} else if (event == 6*K+1){
		++state[0];
	} else { printf("error\n"); }
	return 0; // switch to 1 to stop simulation
}


{

	/* intialize random number generator and seed with timestamp */
	gsl_rng * rng = gsl_rng_alloc (gsl_rng_default);
	#ifdef _RNDSEED
	gsl_rng_set(rng, time(NULL));
	#endif

	/* Gillespie simulation variables */
	int l,i,check;
	double lambda, t, tmp;
	double * rates_data;

    /* Dynamically allocated private arrays must be declared inside 
	 * the parallel region.  */
	#pragma omp parallel shared(rng, inits, record) private(lambda, t, tmp, i, check, rates_data,l)
	{
		/* The vector to store cumulative sum of rates */
		rates_data = (double *) calloc (n_event_types,sizeof(double) );

		/* Loop over ensembles, will be parallelized if compiled with -fopenmp */
		#pragma omp for
		for(l = 0; l < ensembles; l++){
// RESET STATES HERE
			t = 0;
			check=0;
			while(t < max_time){
				/* calculate time until next event */
			rates_calc(rates_data, states, parameters);
			lambda = rates_data[n_event_types-1]; 
				t += gsl_ran_exponential(rng, 1/lambda);
			
				/* Events such as sampling occur at regular intervals */
				fixed_interval_fn(t, states, parameters, record);

				/* Execute the appropriate event */
				tmp = gsl_rng_uniform(rng);
				for(i=0;i<n_event_types;i++){
					if( tmp < rates_data[i]/lambda ){
						check = outcome(state, i);
						break;
					} 
				}
			if(check) break;

		   }/* end evolution */
			free(my_pars);
		}/* end ensembles */
		free(rates_data);
	} /* end parallel */
	gsl_rng_free(rng);
}



