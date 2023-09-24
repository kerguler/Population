#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "population.h"

void custom(const number *qkey, const number num, double *par) {
    par[0] = 480.0 - (qkey[2].i > 4 ? 240.0 : 48.0 * qkey[2].i);
    par[1] = 0.1 * par[0];
}

void fun_transfer(number *key, number num, void *pop) {
    number q[3] = {
        {.i=key[0].i},
        numZERO,
        {.i=key[2].i+1}
    };
    spop2_add(*(population *)pop, q, num);
}

void sim(char stoch) {
    if (stoch)
        spop2_random_init();

    int i, j;

    char arbiters[4] = {AGE_GAMMA, AGE_GAMMA, AGE_CUSTOM, STOP};
    population pop = spop2_init(arbiters, stoch);

    population popdone[3];
    for (i=0; i<3; i++)
        popdone[i] = spop2_init(arbiters, stoch);

    pop->arbiters[0]->fun_q_par = custom;
    pop->numpars[0] = 0;
    //
    pop->arbiters[2]->fun_step = 0;

    number key[3] = {numZERO, numZERO, numZERO};
    number num;
    if (stoch == STOCHASTIC) 
        num.i = 1000;
    else
        num.d = 1000.0;
    spop2_add(pop, key, num);

    if (stoch == STOCHASTIC)
        printf("%d,%d\n",0,0);
    else
        printf("%d,%g\n",0,0.0);

    number size, completed[3];
    double par[2] = {50.0, 10.0};

    for (i=0; i<480; i++) {
        spop2_step(pop, par, &size, completed, popdone);
        if (stoch == STOCHASTIC)
            printf("%d,%d\n",i+1,completed[1].i);
        else
            printf("%d,%g\n",i+1,completed[1].d);

        spop2_foreach(popdone[1], fun_transfer, (void *)(&pop));

        for (j=0; j<3; j++)
            spop2_empty(&popdone[j]);
    }
}

int main(int attr, char *avec[]) {
    if (FALSE)
        sim(DETERMINISTIC);
    else {
        int i;
        for (i=0; i<100; i++)
            sim(STOCHASTIC);
    }

    return 0;
}
