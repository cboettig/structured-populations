/**
 * @file warning_signals.h
 * @author Carl Boettiger, <cboettig@gmail.com>
 * @section DESCRIPTION
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


/*Standard Libraries */
#include <stdio.h>
#include <math.h>
#include <time.h> 
#include <gsl/gsl_randist.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_statistics.h>
#include "gillespie.h"

/* constants provided from main/R */
#define SAMPLE_TIME 50
#define ENSEMBLES 100
#define SAMPLE_FREQ 1.
#define MAX_TIME 500
#define START_POLLUTING 600

/* embedded constants (in intial_conditions() fn)*/
#define No 572 
#define Ao 160
#define N_EVENT_TYPES 2

/** A record of the system state looking back 'sampletime' timesteps*/
typedef struct {
	double * hist;
//	int hist_index;
	double * means;
	double * vars;
	double * skews;
	double * ar1;
	double * arN;
	double * a;
//	int index;

	size_t sampletime;
	double samplefreq;
	size_t windowsize;
	size_t maxtime;
} record;

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


void gslode(void * mypars, double max_time, FILE *theory);
void euler(void *pars, double max_time, FILE *theory);


void warning_signals(
	double * time, 
	double * a, 
	double * means, 
	double * vars,
	double * skews,
	double * ar1, 
	double * arN,	
	int * sample_time, ///< time window over which to compute stats
	double * sample_freq, ///< sampling frequency
	int * max_time, 
	int * n_ensembles, ///< number of ensembles
	double * start_polluting,
	double * pollute_rate,
	double * pollute_increment);

