#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "population.h"

void memorise(const number *qkey, const number num, double *par) {
    par[0] = qkey[1].d + 1.0 / qkey[0].i;
}

void sim( void ) {
    int i;

    char arbiters[3] = {AGE_FIXED, MEMORY, STOP};
    population pop = spop2_init(arbiters, DETERMINISTIC);
    pop->arbiters[1]->fun_q_par = memorise;

    number size, completed;

    number key = numZERO;
    number num = { .d = 1.0 };
    spop2_add(pop, &key, num);

    double par[1] = {10.0};

    for (i=0; i<12; i++) {
        spop2_step(pop, par, &size, &completed, 0);
        spop2_printable(pop,i);
    }

    spop2_free(&pop);
}

int main(int attr, char *avec[]) {
    sim();

    return 0;
}
