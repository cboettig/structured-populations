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
#include <gsl/gsl_statistics.h>
#include "gillespie.h"
#include "pars.h"
#include "record.h"

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

