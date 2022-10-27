#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "population.h"

void run_det(int did, int id, char distr) {
    double par[2] = {40.0, 5.0};

    char arbiters[2] = {distr, STOP};
    population pop = spop2_init(arbiters, DETERMINISTIC);

    number key[1] = {numZERO};
    number num = {.d=10.0};
    spop2_add(pop, key, num);

    printf("%d,%d,%d,%g,%g\n",did,id,0,spop2_size(pop).d,0.0);

    number sz, cm;
    int i;
    for (i=1; i<50; i++) {
        if (i>20) par[0] = 20.0;
        spop2_step(pop, par, &sz, &cm, 0);
        printf("%d,%d,%d,%g,%g\n",did,id,i,sz.d,cm.d);
    }
}

void run_stoch(int did, int id, char distr) {
    double par[2] = {40.0, 5.0};

    char arbiters[2] = {distr, STOP};
    population pop = spop2_init(arbiters, STOCHASTIC);

    number key[1] = {numZERO};
    number num = {.i=10};
    spop2_add(pop, key, num);

    printf("%d,%d,%d,%u,%u\n",did,id,0,spop2_size(pop).i,0);

    number sz, cm;
    int i;
    for (i=1; i<50; i++) {
        if (i>20) par[0] = 20.0;
        spop2_step(pop, par, &sz, &cm, 0);
        printf("%d,%d,%d,%u,%u\n",did,id,i,sz.i,cm.i);
    }
}

int main(int attr, char *avec[]) {
    spop2_random_init();

    char methods[4] = {ACC_ERLANG, ACC_FIXED, AGE_FIXED, AGE_GAMMA};
    char method;

    int i, j;
    for (i=0; i<4; i++) {
        method = methods[i];
        run_det(i, 0, method);
        for (j=0; j<1000; j++)
            run_stoch(i, j+1, method);
    }

    return 0;
}
