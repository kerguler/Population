#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "population.h"

char check(number *key) {
    return key[0].d > 0.5;
}

int main(int attr, char *avec[]) {
    char arbiters[2] = {ACC_ERLANG, STOP};
    population pop = spop2_init(arbiters, DETERMINISTIC);

    double init_key[3] = {0.3, 0.6, 0.9};
    double init_num[3] = {10.0, 15.0, 20.0};
    number key[1];
    number num;

    int i;
    for (i=0; i<3; i++) {
        key[0].d = init_key[i];
        num.d = init_num[i];
        spop2_add(pop, key, num);
    }

    spop2_print(pop);

    printf("Total size: %g\n", spop2_size(pop).d);
    printf("Older than 0.5: %g\n", spop2_count(pop, check).d);

    return 0;
}
