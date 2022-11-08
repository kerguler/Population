#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "population.h"

void run_det(char distr, char mode) {
    double par[2] = {40.0, 5.0};

    char arbiters[2] = {distr, STOP};
    population pop = spop2_init(arbiters, DETERMINISTIC);

    number key[1] = {numZERO};
    number num = {.d=100.0};
    spop2_add(pop, key, num);

    if (mode == 1)
        printf("%d,%g,%g\n",0,spop2_size(pop).d,0.0);
    else if (mode == 2)
        spop2_printable(pop, 0);

    number sz, cm;
    int i;
    for (i=1; i<50; i++) {
        if (i>20) par[0] = 20.0;
        spop2_step(pop, par, &sz, &cm, 0);
        if (mode == 1) 
            printf("%d,%g,%g\n",i,sz.d,cm.d);
        else if (mode == 2)
            spop2_printable(pop, i);
    }
}

int main(int attr, char *avec[]) {
    spop2_random_init();

    if (attr == 2) run_det(ACC_ERLANG, atoi(avec[1]));

    return 0;
}
