#include <gsl/gsl_math.h>
typedef struct {
	double t_step;
	size_t N;
	double * s1;
	double * s2;
} record;

record* record_alloc(size_t N, size_t replicates, double maxtime);
void record_free(record * myrecord);

