#include <stdio.h>
#include <stdlib.h>
#include "population.h"

int main(void) {
    spop2_random_init();

    char arbiters[2] = {ACC_ERLANG, STOP};
    population pop = spop2_init(arbiters, STOCHASTIC);

    number key[1] = {numZERO};
    number num = {.i=10};
    spop2_add(pop, key, num);

    printf("%d,%u,%u\n",0,spop2_size(pop).i,0);

    double par[2] = {20.0, 10.0};

    number sz, cm;
    int i;
    for (i=1; i<50; i++) {
        spop2_step(pop, par, &sz, &cm, 0);
        printf("%d,%u,%u\n",i,sz.i,cm.i);
    }

    return 0;
}