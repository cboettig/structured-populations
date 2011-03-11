#include "optimizers.h"

	typedef struct {
	gsl_vector * vec;
	void * ptr;
    } block;

  block * block_alloc(size_t n)
     {
		block * t = (block *) malloc(sizeof(block));
		t->vec = gsl_vector_alloc(n);
		return t;
     }

  void block_free(block * t)
     {
		gsl_vector_free(t->vec);
		free(t);
     }

	void block_copy(void *inp, void *outp)
	{
		block * in = (block *) inp;
		block * out = (block *) outp;
		gsl_vector_memcpy(out->vec, in->vec);
		out->ptr = in->ptr;
	}

	void * block_copy_construct(void *inp)
	{
		block * x = (block *) inp;
		block * y = block_alloc(x->vec->size);
		block_copy(x, y);
		return y;
	}

	void block_destroy(void *inp)
	{
		block_free( (block *) inp);
	}


	/* distance function */
     double M1(void *xp, void *yp)
     {
       block * x = ((block *) xp);
       block * y = ((block *) yp);
	   int i;
	   double distance = 0;
	   for(i=0; i < x->vec->size; i++){
        	distance += gsl_pow_2(
				gsl_vector_get(x->vec, i) - 
				gsl_vector_get(y->vec, i) 
			);
	   }
	   return sqrt(distance);
     }


	 double round(double);

/* Step function */ 
void S1(const gsl_rng * r, void *xp, double step_size)
{
	block *x = ((block *) xp);
	int i = (int) round(gsl_rng_uniform(r) * (x->vec->size-1));
	gsl_vector_set(x->vec,i,
		//GSL_MAX(0, 
		gsl_vector_get(x->vec,i) + 	gsl_ran_gaussian_ziggurat(r,step_size)
		//)
	); 
}
	 /* print function */
     void P1(void *xp)
     {
		block * myblock = (block *) xp;
		gsl_vector_fprintf(stdout, myblock->vec, "%g");
     }



double E1(void * xp){
	block * myblock = (block *) xp;
	return optim_func(myblock->vec, myblock->ptr);
}



double siman(gsl_vector * vec, void * params, gsl_rng * rng) 
{

	block * myblock = block_alloc(vec->size);
	gsl_vector_memcpy(myblock->vec, vec);
	myblock->ptr = params;


	gsl_siman_params_t siman_params
	   = {N_TRIES, ITERS_FIXED_T, STEP_SIZE,
		  K, T_INITIAL, MU_T, T_MIN};



	gsl_siman_solve(rng, myblock, E1, S1, M1, P1, 
				block_copy, block_copy_construct, 
				block_destroy, 0, siman_params);

	gsl_vector_fprintf(stdout, myblock->vec, "%g");
	printf("siman llik: %g\n", -E1(myblock) );
	P1(myblock);
	return -E1(myblock);
}

