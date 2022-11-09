#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "population.h"

extern gsl_rng *RANDOM;

void sim() {
    char arbiters[3] = {AGE_CONST, AGE_FIXED, STOP};

    population pop = spop2_init(arbiters, DETERMINISTIC);

    {
        number key[2] = {numZERO,numZERO};
        number num = {.d=500.0};
        spop2_add(pop, key, num);
    }

    {
        number key[2] = {{.i=3},{.i=4}};
        number num = {.d=100.0};
        spop2_add(pop, key, num);
    }

    {
        number size, completed[2];
        double par[2] = {0.1, 10};
        spop2_step(pop, par, &size, completed, 0);
    }

    spop2_printable(pop, 0);

    spop2_empty(&pop);

    spop2_printable(pop, 0);
}

int main(int attr, char *avec[]) {
    spop2_random_init();

    sim();

    return 0;
}
