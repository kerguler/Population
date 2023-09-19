#include <stdio.h>
#include <stdlib.h>
#include "population.h"

void run_full() {
    char arbiters[3] = {ACC_ERLANG, ACC_ERLANG, STOP};
    population pop = spop2_init(arbiters, DETERMINISTIC);

    number key[2] = {numZERO,numZERO};
    number num = {.d=10.0};
    spop2_add(pop, key, num);

    printf("%d,%g,%g\n",0,spop2_size(pop).d,0.0);

    double par[4] = {20.0, 10.0, 10.0, 5.0};

    number sz, cm;
    int i;
    for (i=1; i<50; i++) {
        spop2_step(pop, par, &sz, &cm, 0);
        printf("%d,%g,%g\n",i,sz.d,cm.d);
    }
}

void run_save() {
    char arbiters[3] = {ACC_ERLANG, ACC_ERLANG, STOP};
    population pop = spop2_init(arbiters, DETERMINISTIC);

    number key[2] = {numZERO,numZERO};
    number num = {.d=10.0};
    spop2_add(pop, key, num);

    printf("%d,%g,%g\n",0,spop2_size(pop).d,0.0);

    double par[4] = {20.0, 10.0, 10.0, 5.0};

    number sz, cm;
    int i;
    for (i=1; i<20; i++) {
        spop2_step(pop, par, &sz, &cm, 0);
        printf("%d,%g,%g\n",i,sz.d,cm.d);
    }

    printf("\n");
    spop2_printable(pop,0);

    number *state = spop2_savestate(pop);

    population pop2 = spop2_loadstate(state);

    printf("\n");
    spop2_printable(pop2,0);
}

int main(void) {
    run_save();

    return 0;
}