#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <math.h>
#include "population.h"

//Random number variables
const gsl_rng_type *R_TYPE;
gsl_rng *_RAND_GSL;

char *label;

gsl_rng *get_RAND_GSL(void) {
    return _RAND_GSL;
}

void rng_setup(char *lab)
//Can also run like:
//GSL_RNG_TYPE="taus" GSL_RNG_SEED=123 ./GDirect
{
  label = strdup(lab);
  if (_RAND_GSL) {
      fprintf(stderr,"Trying to setup multiple instances of RNG %s\n",label);
      return;
  }
  fprintf(stderr,"Setting up RNG %s\n",label);
  FILE *dev;

  unsigned int i;

  dev=fopen("/dev/urandom","r");

  fread(&i,sizeof(unsigned int),1,dev);

  fclose(dev);

  gsl_rng_env_setup();

  R_TYPE = gsl_rng_mt19937;
  _RAND_GSL = gsl_rng_alloc(R_TYPE);
  gsl_rng_set(_RAND_GSL,i);
}

void rng_setup_seed(unsigned int seed, char *lab)
//Can also run like:
//GSL_RNG_TYPE="taus" GSL_RNG_SEED=123 ./GDirect
{
  label = strdup(lab);
  fprintf(stderr,"Setting up RNG %s with seed = %u\n",label,seed);

  gsl_rng_env_setup();

  R_TYPE = gsl_rng_mt19937;
  _RAND_GSL = gsl_rng_alloc(R_TYPE);
  gsl_rng_set(_RAND_GSL,seed);
}

void rng_destroy()
{
  if (!_RAND_GSL) {
      fprintf(stderr,"RNG %s is already destroyed\n",label);
      return;
  }

  gsl_rng_free(_RAND_GSL);
  _RAND_GSL = 0;

  fprintf(stderr,"RNG %s destroyed\n",label);
}

double rng_exponential(const double mu)
{
  double u = gsl_rng_uniform(_RAND_GSL);
  return -(1.0 / mu) * (u == 0.0 ? 0.0 : log (u));
}
