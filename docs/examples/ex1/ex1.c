#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "population.h"

void run_det(int did, int id, char distr, double *par) {
    char arbiters[2] = {distr, STOP};
    population pop = spop2_init(arbiters, DETERMINISTIC);

    number key[1] = {numZERO};
    number num = {.d=10.0};
    spop2_add(pop, key, num);

    printf("%d,%d,%d,%g,%g\n",did,id,0,spop2_size(pop).d,0.0);

    number sz, cm;
    int i;
    for (i=1; i<50; i++) {
        spop2_step(pop, par, &sz, &cm, 0);
        printf("%d,%d,%d,%g,%g\n",did,id,i,sz.d,cm.d);
    }
}

void run_stoch(int did, int id, char distr, double *par) {
    char arbiters[2] = {distr, STOP};
    population pop = spop2_init(arbiters, STOCHASTIC);

    number key[1] = {numZERO};
    number num = {.i=10};
    spop2_add(pop, key, num);

    printf("%d,%d,%d,%u,%u\n",did,id,0,spop2_size(pop).i,0);

    number sz, cm;
    int i;
    for (i=1; i<50; i++) {
        spop2_step(pop, par, &sz, &cm, 0);
        printf("%d,%d,%d,%u,%u\n",did,id,i,sz.i,cm.i);
    }
}

int main(int attr, char *avec[]) {
    spop2_random_init();

    char methods[7] = {ACC_ERLANG, ACC_PASCAL, ACC_FIXED, AGE_FIXED, AGE_CONST, AGE_GAMMA, AGE_NBINOM};
    char method;

    double par[2];

    int i, j;
    for (i=0; i<7; i++) {
        method = methods[i];
        switch (method) {
            case ACC_ERLANG:
                par[0] = 20.0;
                par[1] = 10.0;
                break;
            case ACC_PASCAL:
                par[0] = 20.0;
                par[1] = 10.0;
                break;
            case ACC_FIXED:
                par[0] = 20.0;
                par[1] = 0.0;
                break;
            case AGE_FIXED:
                par[0] = 20.0;
                par[1] = 0.0;
                break;
            case AGE_CONST:
                par[0] = 1.0/20.0;
                par[1] = 0.0;
                break;
            case AGE_GAMMA:
                par[0] = 20.0;
                par[1] = 10.0;
                break;
            case AGE_NBINOM:
                par[0] = 20.0;
                par[1] = 10.0;
                break;
        }
        run_det(i, 0, method, par);
        for (j=0; j<1000; j++)
            run_stoch(i, j+1, method, par);
    }

    return 0;
}
