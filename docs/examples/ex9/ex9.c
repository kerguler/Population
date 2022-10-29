#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "population.h"

extern gsl_rng *RANDOM;

void sim(char stoch) {
    int i, j;

    double death[5] = {0.07, 0.004, 0.003, 0.0025, 0.27};
    double develop[5] = {0.6, 5.0, 5.9, 4.1, 1e13};
    double A0 = 600.0;
    double q = 8.5;
    double tau = 24.0;

    char arbiters[3] = {AGE_CONST, AGE_FIXED, STOP};

    population pop[5];
    for (i=0; i<5; i++)
        pop[i] = spop2_init(arbiters, stoch);

    population popdone[5][2];
    for (i=0; i<5; i++)
        for (j=0; j<2; j++)
            popdone[i][j] = spop2_init(arbiters, stoch);

    number key[2] = {numZERO,numZERO};
    number num;
    if (stoch == STOCHASTIC) 
        num.i = 500;
    else
        num.d = 500.0;
    spop2_add(pop[0], key, num);

    if (stoch == STOCHASTIC) {
        printf("%d",0); for (i=0; i<5; i++) printf(",%d",spop2_size(pop[i]).i); printf("\n");
    } else {
        printf("%d",0); for (i=0; i<5; i++) printf(",%g",spop2_size(pop[i]).d); printf("\n");
    }

    number size[5], completed[5][2];
    double par[2];
    member elm, tmp;
    number eggs;

    for (j=0; j<300*tau; j++) {
        for (i=0; i<5; i++) {
            par[0] = death[i] / tau;
            par[1] = develop[i] * tau;
            spop2_step(pop[i], par, &size[i], completed[i], popdone[i]);

            if (i) {
                HASH_ITER(hh, popdone[i-1][1]->members, elm, tmp) {
                    spop2_add(pop[i], key, elm->num);
                }
            }
        }

        if (stoch == STOCHASTIC) {
            eggs.i = q * size[4].i * exp(-size[4].i / A0) / tau;
            eggs.i = gsl_ran_poisson(RANDOM, eggs.i);
        } else {
            eggs.d = q * size[4].d * exp(-size[4].d / A0) / tau;
        }
        spop2_add(pop[0], key, eggs);

        if (stoch == STOCHASTIC) {
            printf("%d",0); for (i=0; i<5; i++) printf(",%d",spop2_size(pop[i]).i); printf("\n");
        } else {
            printf("%d",0); for (i=0; i<5; i++) printf(",%g",spop2_size(pop[i]).d); printf("\n");
        }
    }
}

int main(int attr, char *avec[]) {
    spop2_random_init();

    if (TRUE)
        sim(DETERMINISTIC);
    else {
        int i;
        for (i=0; i<100; i++)
            sim(STOCHASTIC);
    }

    return 0;
}
