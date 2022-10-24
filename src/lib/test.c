#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "population.h"

void test() {
    char arbiters[3] = {AGE_GAMMA, ACC_ERLANG, STOP};
    population pop = spop2_init(arbiters, DETERMINISTIC);

    number key[2] = {{.i=10}, {.d=0.2}};
    number num = {.d=10.0};
    spop2_add(pop, key, num);

    number keyb[2] = {{.i=5}, {.d=0.2}};
    number numb = {.d=6.0};
    spop2_add(pop, keyb, numb);

    number key2[2] = {{.i=20}, {.d=0.01}};
    number num2 = {.d=5.0};
    spop2_add(pop, key2, num2);

    spop2_print(pop);

    double h, hz;
    hazpar hp = pop->arbiters[0]->fun_pars(10.0, 2.0);
    printf("\nk = %g\ntheta = %g\nstay = %d\n", hp.k.d, hp.theta, hp.stay);

    unsigned int i = 8;
    h = pop->arbiters[0]->fun_haz(i, hp.k, hp.theta);
    printf("h = %g\n", h);

    number q = {.i=7};
    hz = pop->arbiters[0]->fun_calc(pop->arbiters[0]->fun_haz, 1, q, hp.k, hp.theta, key);
    printf("hz = %g\n", hz);

    number completed[2];
    double par[4] = {10.0, 2.0, 9.0, 3.0};
    spop2_step(pop, par, completed, 0);

    spop2_print(pop);
    printf("Completed: %g %g\n",completed[0].d,completed[1].d);

    spop2_free(&pop);

    number a = {.d=random()};
    printf("%g %u\n",a.d,a.i);
    a.i = 0;
    printf("%g %u\n",a.d,a.i);
}

void test2() {
    char arbiters[2] = {ACC_ERLANG, STOP};
    population pop = spop2_init(arbiters, STOCHASTIC);

    number key[1] = {{.d=0}};
    number num = {.i=10};
    spop2_add(pop, key, num);

    printf("%d %u\n",0,spop2_size(pop).i);

    double par[2] = {10.0, 2.0};
    int i;
    for (i=1; i<20; i++) {
        spop2_step(pop, par, 0, 0);
        printf("%d %u\n",i,spop2_size(pop).i);
    }

    spop2_print(pop);
}

int main(int attr, char *avec[]) {
    random_init();

    test2();

    return 0;
}
