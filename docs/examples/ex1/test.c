#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "population.h"

void run_det(int id, char distr, double *par) {
    char arbiters[2] = {distr, STOP};
    population pop = spop2_init(arbiters, DETERMINISTIC);

    number key[1] = {numZERO};
    number num = {.d=10.0};
    spop2_add(pop, key, num);

    printf("%d,%d,%g,%g\n",id,0,spop2_size(pop).d,0.0);

    number sz, cm;
    int i;
    for (i=1; i<50; i++) {
        spop2_step(pop, par, &sz, &cm, 0);
        printf("%d,%d,%g,%g\n",id,i,sz.d,cm.d);
    }
}

void run_stoch(int id, char distr, double *par) {
    char arbiters[2] = {distr, STOP};
    population pop = spop2_init(arbiters, STOCHASTIC);

    number key[1] = {numZERO};
    number num = {.i=10};
    spop2_add(pop, key, num);

    printf("%d,%d,%u,%u\n",id,0,spop2_size(pop).i,0);

    number sz, cm;
    int i;
    for (i=1; i<50; i++) {
        spop2_step(pop, par, &sz, &cm, 0);
        printf("%d,%d,%u,%u\n",id,i,sz.i,cm.i);
    }
}

int main(int attr, char *avec[]) {
    random_init();

    char method = AGE_NBINOM;
    double par[2] = {20.0, 10.0};
    //double par[1] = {1.0/20.0};
    //double par[1] = {20.0};

    run_det(0, method, par);
    int i;
    for (i=0; i<1000; i++)
        run_stoch(i+1, method, par);

    return 0;
}
