#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "population.h"

double custom(hazard hfun, unsigned int d, number q, number k, double theta, const number *qkey) {
    double devmn = 480.0 - (qkey[2].i > 4 ? 240.0 : 48.0 * qkey[2].i);
    double devsd = 0.1 * devmn;
    hazpar hz = age_gamma_pars(devmn, devsd);
    return age_hazard_calc(0, 0, qkey[0].i, hz.k, hz.theta, qkey);
}

void sim() {
    char arbiters[2] = {AGE_CUSTOM, AGE_GAMMA, AGE_DUMMY, STOP};
    population pop = spop2_init(arbiters, DETERMINISTIC);

    pop->arbiters[0]->fun_calc = custom;

    number key[1] = {numZERO};
    number num = {.d=100.0};
    spop2_add(pop, key, num);

    printf("%d,%g,%g\n",0,spop2_size(pop).d,0.0);

    number size, completed;

    int i;
    for (i=0; i<100; i++) {
        spop2_step(pop, 0, &size, &completed, 0);
        printf("%d,%g,%g\n",i+1,size.d,completed.d);
    }
}

int main(int attr, char *avec[]) {
    spop2_random_init();

    sim();

    return 0;
}
