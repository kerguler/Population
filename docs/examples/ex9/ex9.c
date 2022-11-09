#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "population.h"

extern gsl_rng *RANDOM;

void sim(char stoch, char id) {
    int i, j, k;

    double death[5] = {0.07, 0.004, 0.003, 0.0025, 0.27};
    double develop[5] = {0.6, 5.0, 5.9, 4.1, 1e13};
    double A0 = 600.0;
    double q = 8.5;
    double tau = 24.0;

    char arbiters[3] = {AGE_CONST, AGE_FIXED, STOP};

    population pop[5];
    for (i=0; i<5; i++)
        pop[i] = spop2_init(arbiters, stoch);

    number key[2] = {numZERO,numZERO};
    number num;
    if (stoch == STOCHASTIC) 
        num.i = 500;
    else
        num.d = 500.0;
    spop2_add(pop[0], key, num);

    if (stoch == STOCHASTIC) {
        printf("%d,%d",id,0); for (i=0; i<5; i++) printf(",%d",spop2_size(pop[i]).i); printf("\n");
    } else {
        printf("%d,%d",id,0); for (i=0; i<5; i++) printf(",%g",spop2_size(pop[i]).d); printf("\n");
    }

    number size[5], completed[5][2];
    double par[2];
    number eggs;

    for (j=0; j<300*tau; j++) {
        for (i=0; i<5; i++) {
            par[0] = death[i] / tau;
            par[1] = develop[i] * tau;
            spop2_step(pop[i], par, &size[i], completed[i], 0);

            if (i) {
                spop2_add(pop[i], key,completed[i-1][1]);
            }
        }

        if (stoch == STOCHASTIC) {
            double tmp = q * (double)size[4].i * exp(-(double)size[4].i / A0) / tau;
            eggs.i = (unsigned int)gsl_ran_poisson(RANDOM, tmp);
        } else {
            eggs.d = q * size[4].d * exp(-size[4].d / A0) / tau;
        }
        spop2_add(pop[0], key, eggs);

        if (stoch == STOCHASTIC) {
            printf("%d,%d",id,j); for (k=0; k<5; k++) printf(",%d",spop2_size(pop[k]).i); printf("\n");
        } else {
            printf("%d,%d",id,j); for (k=0; k<5; k++) printf(",%g",spop2_size(pop[k]).d); printf("\n");
        }
    }
}

int main(int attr, char *avec[]) {
    spop2_random_init();

    if (FALSE)
        sim(DETERMINISTIC,0);
    else {
        int i;
        for (i=0; i<100; i++)
            sim(STOCHASTIC,i+1);
    }

    return 0;
}
