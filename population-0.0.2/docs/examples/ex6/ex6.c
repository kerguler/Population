#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "population.h"

void sim() {
    char arbiters[3] = {AGE_CONST, ACC_ERLANG, STOP};
    population pop = spop2_init(arbiters, DETERMINISTIC);

    number key[2] = {numZERO,numZERO};
    number num = {.d=100.0};
    spop2_add(pop, key, num);

    printf("0,%d,%g,%g,%g\n",0,spop2_size(pop).d,0.0,0.0);

    number size, completed[2];
    population popdone[2];
    popdone[0] = spop2_init(arbiters, DETERMINISTIC);
    popdone[1] = spop2_init(arbiters, DETERMINISTIC);
    double par[4] = {1.0/60.0, 30.0, 5.0};

    int i;
    for (i=0; i<100; i++) {
        spop2_step(pop, par, &size, completed, popdone);
        printf("0,%d,%g,%g,%g\n",i+1,size.d,completed[0].d,completed[1].d);
        spop2_printable(popdone[0], 1);
        spop2_printable(popdone[1], 2);
        spop2_empty(&popdone[0]);
        spop2_empty(&popdone[1]);
    }
}

int main(int attr, char *avec[]) {
    spop2_random_init();

    sim();

    return 0;
}
