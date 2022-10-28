#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "population.h"

double custom(hazard hfun, unsigned int d, number q, number k, double theta, const number *qkey) {
    return 1.0/20.0;
}

void sim() {
    char arbiters[2] = {AGE_CUSTOM, STOP};
    population pop = spop2_init(arbiters, DETERMINISTIC);

    pop->arbiters[0]->fun_calc = custom;

    number key[1] = {numZERO};
    number num = {.d=100.0};
    spop2_add(pop, key, num);

    printf("%d,%g,%g\n",0,spop2_size(pop).d,0.0);

    number size, completed;
    double par[2] = {0.0, 0.0};

    int i;
    for (i=0; i<100; i++) {
        spop2_step(pop, par, &size, &completed, 0);
        printf("%d,%g,%g\n",i+1,size.d,completed.d);
    }
}

int main(int attr, char *avec[]) {
    spop2_random_init();

    sim();

    return 0;
}
