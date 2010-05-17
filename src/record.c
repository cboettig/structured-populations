#include "record.h"

record * record_alloc(size_t sampletime, double samplefreq, size_t maxtime)
{
	size_t windowsize = sampletime/samplefreq;
	size_t length = maxtime/samplefreq;
	record * my_record = (record *) malloc(sizeof(record) );
	my_record->hist  = (double *) calloc( windowsize, sizeof(double) );
	my_record->means = (double *) calloc( length, sizeof(double) );
	my_record->vars  = (double *) calloc( length, sizeof(double) );
	my_record->skews = (double *) calloc( length, sizeof(double) );
	my_record->ar1   = (double *) calloc( length, sizeof(double) );
	my_record->arN   = (double *) calloc( length, sizeof(double) );
//	my_record->index = 0;
	my_record->windowsize = windowsize;
	my_record->sampletime = sampletime;
	my_record->samplefreq = samplefreq;
	my_record->maxtime = maxtime;
//	my_record->hist_index = 0;
	my_record->a   = (double *) calloc( length, sizeof(double) );
	return my_record;
}

void record_free(record * my_record)
{
	free(my_record->hist);
	free(my_record->means);
	free(my_record->vars);
	free(my_record->skews);
	free(my_record->ar1);
	free(my_record->arN);
	free(my_record->a);
	free(my_record);
}
