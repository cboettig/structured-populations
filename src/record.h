#include <stdlib.h>

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

record * record_alloc(size_t sampletime, double samplefreq, size_t maxtime);
void record_free(record * my_record);

