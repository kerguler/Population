#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "population.h"

double fun_harvest(number *key) {
    return key[0].i > 136 ? 1.0 : 0.0;
}

double fun_rest(number *key) {
    return 1.0;
}

void sim(char stoch) {
    int i;

    char arbiters[2] = {AGE_GAMMA, STOP};
    population pop[3];
    pop[0] = spop2_init(arbiters, stoch);
    pop[1] = spop2_init(arbiters, stoch);
    pop[2] = spop2_init(arbiters, stoch);
    population popdone = spop2_init(arbiters, stoch);

    number key = numZERO;
    number num;
    if (stoch == STOCHASTIC) 
        num.i = 10;
    else
        num.d = 10.0;
    spop2_add(pop[0], &key, num);

    if (stoch == STOCHASTIC)
        printf("%d,%d,%d,%d\n",0,0,0,0);
    else
        printf("%d,%g,%g,%g\n",0,0.0,0.0,0.0);

    number size, completed;
    double par[2] = {120.0, 24.0};

    for (i=0; i<240; i++) {
        spop2_step(pop[0], par, &size, &completed, &popdone);
        spop2_harvest(popdone, pop[1], fun_harvest);
        spop2_harvest(popdone, pop[2], fun_rest);

        if (stoch == STOCHASTIC)
            printf("%d,%d,%d,%d\n",i+1,spop2_size(pop[0]).i,spop2_size(pop[1]).i,spop2_size(pop[2]).i);
        else
            printf("%d,%g,%g,%g\n",i+1,spop2_size(pop[0]).d,spop2_size(pop[1]).d,spop2_size(pop[2]).d);

        spop2_empty(&popdone);
    }
}

int main(int attr, char *avec[]) {
    if (FALSE)
        sim(DETERMINISTIC);
    else {
        spop2_random_init();
        int i;
        for (i=0; i<100; i++)
            sim(STOCHASTIC);
    }

    return 0;
}
