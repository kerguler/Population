#include "population.h"

#define NumPar 4
#define NumMet 2
#define NumEnv 1

extern gsl_rng *RANDOM;

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
    "female", "egg",
    "initf", "inite", "mean", "std"
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

void sim(double               *envar,
         double               *param,
         int                  *ftime,
         int                 *repeat,
         double              *result,
         int                *success) {

    char arbiters[3] = {ACC_ERLANG, ACC_ERLANG, STOP};
    population adult = spop2_init(arbiters, STOCHASTIC);
    population poptable[2] = {
        spop2_init(arbiters, STOCHASTIC),
        spop2_init(arbiters, STOCHASTIC)
    };

    double par[4] = {param[0], param[1], param[2], param[3]};
    double fec = param[4];
    number size = numZERO;
    number completed[2] = {numZERO, numZERO};

    unsigned int eggs = 0;
    int rm = 0, tm = 0;

    for (rm=0; rm < *repeat; rm++) {
        eggs = 0;
        tm = 0;
        print_out(rm, tm, *ftime, spop2_size(adult).i, eggs, result);
        for (tm=1; tm < *ftime; tm++) {
            spop2_step(adult, par, &size, completed, poptable);

            eggs = gsl_ran_poisson(RANDOM, completed[1].i * fec);

            spop2_addpop(adult, poptable[1]);

            print_out(rm, tm, *ftime, size.i, eggs, result);
        }
        spop2_empty(&adult);
        spop2_empty(&poptable[0]);
        spop2_empty(&poptable[1]);
    }

    spop2_free(&adult);
    spop2_free(&poptable[0]);
    spop2_free(&poptable[1]);
}
