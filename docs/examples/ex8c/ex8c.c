#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "population.h"

#define N 5
#define NEWAGE TRUE
#define STOCH FALSE

void fun_harvest(number *key, number num, number *newkey, double *frac) {
    newkey[0].i = NEWAGE ? 0 : key[0].i;
    *frac = 0.5;
}

void fun_rest(number *key, number num, number *newkey, double *frac) {
    newkey[0].i = key[0].i;
    *frac = 1.0;
}

void sim(char stoch) {
    int i, j;

    char arbiters[2] = {AGE_GAMMA, STOP};
    population *pop = (population *)malloc(N * sizeof(population));
    population *popdone = (population *)malloc(N * sizeof(population));
    for (j=0; j<N; j++) {
        pop[j] = spop2_init(arbiters, stoch);
        popdone[j] = spop2_init(arbiters, stoch);
    }

    number key = numZERO;
    number num;
    if (stoch == STOCHASTIC) 
        num.i = 100;
    else
        num.d = 100.0;
    spop2_add(pop[0], &key, num);

    printf("0");
    if (stoch == STOCHASTIC)
        for (j=0; j<N; j++)
            printf(",0");
    else
        for (j=0; j<N; j++)
            printf(",0.0");
    printf("\n");

    number size[N], completed[N];
    double par[2] = {48.0, 6.0};

    for (i=0; i<240; i++) {
        for (j=0; j<(N-2); j++) {
            spop2_step(pop[j], par, &size[j], &completed[j], &popdone[j]);
        }
        for (j=0; j<(N-2); j++) {
            spop2_harvest(popdone[j], pop[j+1], fun_harvest);
            spop2_harvest(popdone[j], pop[j+2], fun_rest);
        }

        printf("%d",i+1);
        if (stoch == STOCHASTIC)
            for (j=0; j<N; j++)
                printf(",%d",spop2_size(pop[j]).i);
        else
            for (j=0; j<N; j++)
                printf(",%g",spop2_size(pop[j]).d);
        printf("\n");

        for (j=0; j<N; j++) {
            spop2_empty(&popdone[j]);
        }
    }

    for (j=0; j<N; j++) {
        spop2_free(&pop[j]);
        spop2_free(&popdone[j]);
    }
}

int main(int attr, char *avec[]) {
    if (STOCH) {
        spop2_random_init();
        int i;
        for (i=0; i<100; i++)
            sim(STOCHASTIC);
    } else
        sim(DETERMINISTIC);

    return 0;
}
