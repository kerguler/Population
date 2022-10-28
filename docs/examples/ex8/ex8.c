#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "population.h"

double custom(hazard hfun, unsigned int d, number q, number k, double theta, const number *qkey) {
    double devmn = 480.0 - (qkey[2].i > 4 ? 240.0 : 48.0 * qkey[2].i);
    double devsd = 0.1 * devmn;
    hazpar hz = age_gamma_pars(devmn, devsd);
    double a = age_hazard_calc(age_gamma_haz, 0, qkey[0], hz.k, hz.theta, qkey);
    return a;
}

void sim_det(char stoch) {
    int i, j;

    char arbiters[4] = {AGE_CUSTOM, AGE_GAMMA, AGE_CUSTOM, STOP};
    population pop = spop2_init(arbiters, stoch);

    population popdone[3];
    for (i=0; i<3; i++)
        popdone[i] = spop2_init(arbiters, stoch);

    pop->arbiters[0]->fun_calc = custom;
    pop->arbiters[2]->fun_step = 0;

    number key[3] = {numZERO,numZERO,numZERO};
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

    if (FALSE)
        sim_det(DETERMINISTIC);
    else {
        int i;
        for (i=0; i<100; i++)
            sim_det(STOCHASTIC);
    }

    return 0;
}
