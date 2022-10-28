#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "population.h"

void sim(char stoch) {
    int i, j;

    double death[5] = {0.07, 0.004, 0.003, 0.0025, 0.27};
    double develop[5] = {0.6, 5.0, 5.9, 4.1, 1e13};
    double A0 = 600.0;
    double q = 8.5;
    double tau = 24.0;

    char arbiters[3] = {AGE_CONST, AGE_FIXED, STOP};

    population pop[5];
    for (i=0; i<5; i++)
        pop[i] = spop2_init(arbiters, stoch);

    population popdone[2];
    for (i=0; i<2; i++)
        popdone[i] = spop2_init(arbiters, stoch);

    number key[2] = {numZERO,numZERO};
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

    member elm, tmp;

    number q[3] = {numZERO,numZERO,numZERO};

    for (i=0; i<480; i++) {
        spop2_step(pop, par, &size, completed, popdone);
        if (stoch == STOCHASTIC)
            printf("%d,%d\n",i+1,completed[1].i);
        else
            printf("%d,%g\n",i+1,completed[1].d);
        //
        HASH_ITER(hh, popdone[1]->members, elm, tmp) {
            q[0].i = elm->key[0].i+1;
            q[1] = numZERO;
            q[2].i = elm->key[2].i+1;
            spop2_add(pop, q, elm->num);
            for (j=0; j<3; j++)
                spop2_empty(&popdone[j]);
        }
    }
}

int main(int attr, char *avec[]) {
    spop2_random_init();

    if (TRUE)
        sim_det(DETERMINISTIC);
    else {
        int i;
        for (i=0; i<100; i++)
            sim_det(STOCHASTIC);
    }

    return 0;
}
