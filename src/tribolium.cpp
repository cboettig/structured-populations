/**
 * @file tribolium.cpp
 * @author Carl Boettiger, <cboettig@gmail.com>
 * @section DESCRIPTION
 * This is the actual population dynamics simulation.  Currently 
 * completely contained in this file except for the dependence on
 * my linked list library.  Implements an individual-based, Gillespie
 * simulation which captures demographic stochasticity in continuous 
 * time through four different age classes (egg, larvae, pupae, adult).
 * 
 * @section LICENCE
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License as
 * published by the Free Software Foundation; either version 3 of
 * the License, or (at your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
 * General Public License for more details at
 * http://www.gnu.org/copyleft/gpl.html
 *
 */


#include <stdio.h>
#include <math.h>
#include <time.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include "kde.h"

/*
 * Events:
 * -- random times --
 * death of an egg
 * death of larva
 * death of pupa
 * death of adult
 * birth of egg
 *
 *  -- fixed times --
 * transition egg to larva
 * transition larva to pupa
 * transition pupa to adult
 *
 * compare time to next event to time to next maturation
 */

#include <iostream>
#include <list>
using namespace std;
extern "C" {

/** Simulate function, not directly callable from R */
void simulate(int *state, double *pars, double *dt, double *T, gsl_rng *rng, FILE * file)
{

	list<double> pop;
	list<double>::iterator p;

	/* Biological Parameters.  eggs and pupa can get canibalized */
	double b = pars[0];											/* Birth rate (per day) */
	double u_egg = pars[1], u_larva = pars[2], u_pupa = pars[3], u_adult = pars[4]; /* Mortality rate (per day) */
	double a_egg = pars[5], a_larva = pars[6], a_pupa = pars[7];	/* age at which each stage matures, in days */
	double cannibal_larva_eggs = pars[8], cannibal_adults_pupa = pars[9], cannibal_adults_eggs = pars[10]; /* cannibalism of x on y, per day */
	double a_larva_asym = pars[11];									/* age after which larval size asymptotes */


	/* Gillespie simulation parameters */
	double t=0, sampletime = *dt, timestep = 0;					/* overall time, sample timer, random time increment */
	int event=-1, birthdeath;
	double tmp;
	double rates[5], cumrates[5], sumrates = 0;
	int age_classes[5];

	/* some additional parameters we keep track of */
	double immature_larva_weights;
	double next_mature_time;


	/* initialize No individuals */
	int i,j;
	for(i=0;i<state[0];i++) pop.push_back(0.0);
	for(i=0;i<state[1];i++)	pop.push_back(a_egg);
	for(i=0;i<state[2];i++)	pop.push_back(a_larva_asym);
	for(i=0;i<state[3];i++) pop.push_back(a_larva);
	for(i=0;i<state[4];i++)	pop.push_back( a_pupa);
	
	/* initialize simulation */
	t = 0, timestep = 0, sumrates = 0;
	while(t <= *T){
		for(j=0;j<5;j++) age_classes[j] = 0;
		immature_larva_weights = 0;
		next_mature_time = a_pupa+100;
		/** Tour the population.  Keep track of:
		 * (a) number falling into each age bracket 
		 * (b) if any will mature before the next birth/death event
		 * (c) get a pointer to a member of each class, so we can grab
		 *     a member of the appropriate class when we have to kill it 
		 *
		 *	immature larva have a weight determined by their age, so we'll calculate that too
		 *     */
		for(p = pop.begin(); p != pop.end(); p++) {
			*p += timestep;
			if(*p < a_egg)
			{
				++age_classes[0];												
				next_mature_time = GSL_MIN_DBL(next_mature_time, a_egg - *p);
			} else if(*p < a_larva_asym) {
				++age_classes[1];
				immature_larva_weights += (*p - a_egg)/(a_larva_asym-a_egg);
				next_mature_time = GSL_MIN_DBL(next_mature_time, a_larva_asym - *p);
			} else if(*p < a_larva) {
				++age_classes[2]; 
				next_mature_time = GSL_MIN_DBL(next_mature_time, a_larva - *p);
			} else if(*p < a_pupa) { 
				++age_classes[3]; 
				next_mature_time = GSL_MIN_DBL(next_mature_time, a_pupa - *p);
			} else if (*p >= a_pupa) {  
				++age_classes[4];
			}
		}
		/* Create the Gillespie vector of rates, in this order: 
		 * egg death, larva death, pupa death, adult death, birth 
		 * Then determine which event happens and when */
		rates[0] = age_classes[0] * (u_egg + ((double) age_classes[2] + immature_larva_weights) * cannibal_larva_eggs + age_classes[4]*cannibal_adults_eggs );
		rates[1] = (age_classes[2]+age_classes[1])*u_larva;
		rates[2] = age_classes[3]* (u_pupa + age_classes[4]*cannibal_adults_pupa );
		rates[3] = age_classes[4]*u_adult;
		rates[4] = age_classes[4]*b;

		sumrates=0;
		for(i=0;i<5;i++){
			sumrates += rates[i];
			cumrates[i] = sumrates;
		}

		timestep = gsl_ran_exponential(rng, 1/sumrates);
		birthdeath = 0;
		/* If no one is maturing before next event, pick a random birth/death */
		if(timestep < next_mature_time){
			birthdeath = 1;
		} else { /* next event is someone maturing */
			timestep = next_mature_time;
		}
		t += timestep;
		/* Sample the data at regular intervals */
		while(t > sampletime && sampletime < *T){
			/* update state */
			for(i=0;i<5;i++) state[i] = age_classes[i];
			fprintf(file, "%lf %d %d %d %d %d\n", sampletime, state[0], state[1], state[2], state[3], state[4]);
			sampletime += *dt;
		}

		if(birthdeath){
			tmp = gsl_rng_uniform(rng);
			for(i=0;i<5;i++){
				if (tmp < cumrates[i]/sumrates){
					event = i;
					break;
				}
			}
			for(p = pop.begin(); p != pop.end(); p++) {
				if(*p < a_egg && event == 0){ pop.erase(p); break;
				} else if(*p < a_larva_asym && event == 1){ pop.erase(p); break;
				} else if(*p < a_larva && event == 1){ pop.erase(p); break;
				} else if(*p < a_pupa && event == 2){  pop.erase(p); break;
				} else if (*p >= a_pupa && event == 3){  pop.erase(p); break; 
				}
			}
			if(event == 4) pop.push_back(0.0);
		}
	} // loop over time
	for(i=0;i<5;i++) state[i] = age_classes[i];
	pop.clear();
}


/** R wrapper to simulate */
void beetle_sim(int *state, double *pars, double *dt, double *T, int *seed){
	gsl_rng *rng = gsl_rng_alloc( gsl_rng_default);
	gsl_rng_set(rng, *seed);
	FILE * file = fopen("beetle_sim.txt", "w");
	simulate(state, pars, dt, T, rng, file);
	gsl_rng_free(rng);
	fclose(file);
}


/* Simulate an ensemble in order to estimate probability of an observation, can be called directly from R*/
void ensemble(int *state, int *initial, double *pars, double *dt, int *seed, int *reps, double *probs, int *nstates)
{
	/* Initialize random number generator */
	gsl_rng *rng = gsl_rng_alloc( gsl_rng_default);
	gsl_rng_set(rng, *seed);

	int i,j;
	/* preserve the original initial conditions */
	int *orig = (int *) malloc(*nstates*sizeof(int));
	for(j = 0; j < *nstates; j++)
	{
		orig[j] = initial[j];
	}

	/* perform the simulations */
	double *replicates = (double *) calloc(*nstates * *reps, sizeof(double));
	for(i = 0; i < *reps; i++)
	{
		simulate(initial, pars, dt, dt, rng, stdout);
		for(j = 0; j < *nstates; j++)
		{
			replicates[*reps * j + i] = (double) initial[j];
			initial[j] = orig[j];
		}
	}

	/* estimate the kernel density */
	for(j = 0; j < *nstates; j++)
	{
		probs[j] = kerneldensity(&replicates[*reps * j], (double) state[j], *reps);
	}
	free(orig);
	free(replicates);
	gsl_rng_free(rng);
}

}



extern "C" {
int tibol(void){
	const double b = 5;		// Birth rate (per day) 
	const double u_egg = 0.00001, 
				 u_larva = .001,
				 u_pupa = 0.0,
				 u_adult = 0.003; // Mortality rate (per day) 
	const double a_egg = 3.8,
				 a_larva = 20.2,
				 a_pupa = 25.2;	// age at which each stage matures, in days 
	const double cannibal_larva_eggs = 0.01,
				 cannibal_adults_pupa = 0.004,
				 cannibal_adults_eggs =  0.01; // cannibalism of x on y, per day
	const double a_larva_asym = 11.8;	  // age after which larval size asymptotes 
	double pars[12] = {	b, u_egg, u_larva, u_pupa, u_adult, a_egg,
						a_larva, a_pupa, cannibal_larva_eggs,
						cannibal_adults_pupa, cannibal_adults_eggs,
						a_larva_asym};
	int nstates = 5;
	int initial[] = {450, 0, 0, 0, 30};
	int state[] = {450, 0, 0, 1, 31};
	int seed = time(NULL);
	double dt = 14;
	int reps = 20;
	double probs[5];
	double T = 200;

	ensemble(state, initial, pars, &dt, &seed, &reps, probs, &nstates);
	printf("%g %g %g %g %g\n", probs[0], probs[1], probs[2], probs[3], probs[4]);

	beetle_sim(initial, pars, &dt, &T, &seed); 

	return 0;
}

} // ends the  extern "C"


