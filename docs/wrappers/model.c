#include "population.h"

extern gsl_rng *RANDOM;

#define NumPar 7
#define NumMet 2
#define NumEnv 1

void init(int *np, int *nm) {
  spop2_random_init();

  *np = NumPar;
  *nm = NumMet;
}

void destroy(void) {
  spop2_random_destroy();  
}

void parnames(char **names) {
  char temp[NumMet+NumPar][256] = {
    "adult", "egg",
    "init.adult", "init.egg", "lifetime.mn", "lifetime.sd", "gon.cycle.mn", "gon.cycle.sd", "fecundity"
  };
  int i;
  for (i=0; i<(NumMet+NumPar); i++)
    names[i] = strdup(temp[i]);
}

void print_out(int rep, int tm, int ftime,
               unsigned int female, unsigned int egg,
               double *result) {
    result += (rep * ftime * NumMet) + tm * NumMet;
    result[0] = (double)female;
    result[1] = (double)egg;
}

void fun_transfer(number *key, number num, void *pop) {
    number elmkey[2] = {key[0], numZERO};
    spop2_add(*((population *)pop), elmkey, num);
}

void sim(double               *envar,
         double               *param,
         int                  *ftime,
         int                 *repeat,
         double              *result,
         int                *success) {

    *success = 0;

    char arbiters[3] = {ACC_ERLANG, ACC_ERLANG, STOP};
    population adult = spop2_init(arbiters, STOCHASTIC);
    population poptable[2] = {
        spop2_init(arbiters, STOCHASTIC),
        spop2_init(arbiters, STOCHASTIC)
    };

    number init[2] = {{.i=param[0]},{.i=param[1]}};
    double par[4] = {param[2], param[3], param[4], param[5]};
    double fec = param[6];
    number size = numZERO;
    number key[2] = {numZERO, numZERO};
    number completed[2] = {numZERO, numZERO};

    unsigned int eggs = 0;
    int rm = 0, tm = 0;

    for (rm=0; rm < *repeat; rm++) {
        tm = 0;
        spop2_add(adult, key, init[0]);
        eggs = init[1].i;
        print_out(rm, tm, *ftime, spop2_size(adult).i, eggs, result);

        for (tm=1; tm < *ftime; tm++) {
            spop2_step(adult, par, &size, completed, poptable);

            eggs = gsl_ran_poisson(RANDOM, completed[1].i * fec);

            spop2_foreach(poptable[1], fun_transfer, (void *)(&adult));

            print_out(rm, tm, *ftime, spop2_size(adult).i, eggs, result);

            spop2_empty(&poptable[0]);
            spop2_empty(&poptable[1]);
        }

        spop2_empty(&adult);
    }

    spop2_free(&adult);
    spop2_free(&poptable[0]);
    spop2_free(&poptable[1]);

    *success = 1;
}
