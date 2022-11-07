#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "population.h"

extern gsl_rng *RANDOM;

#define briere1(T,T1,T2,a) ((a)*(T)*((T)-(T1))*sqrt((T2)-(T)))

double var_temp(double mu, double sd) {
    if (sd == 0.0) return mu;
    double y = mu + gsl_ran_gaussian(RANDOM,sd);
    y = min(50.0, max(0.0, y));
    return y;
}

double dev_time(double temp) {
    double y = briere1(temp, 0, 50, 1.5e-5);
    return y == 0.0 ? 1e3 : 1.0 / y;
}

void sim(double mu, double sd, int rep, double eps, char verbose) {
    spop2_set_eps(eps);
    //
    double temp, dev;
    number size, completed;
    double par[2] = {0.0, 5.0};
    char arbiters[2] = {ACC_ERLANG, STOP};
    number key[1] = {numZERO};
    number num = {.d=100.0};
    population pop;
    int r, i;
    //
    for (r=0; r<rep; r++) {
        if (verbose == 1) printf("Rep %d\n",r);

        pop = spop2_init(arbiters, DETERMINISTIC);

        spop2_add(pop, key, num);

        printf("%g,%g,%d,%d,NaN,%g,%g\n",mu,eps,r,0,spop2_size(pop).d,0.0);

        for (i=0; i<100; i++) {
            temp = var_temp(mu, sd);
            dev = dev_time(temp);

            par[0] = dev;
            spop2_step(pop, par, &size, &completed, 0);

            printf("%g,%g,%d,%d,%g,%g,%g\n",mu,eps,r,i+1,dev,size.d,completed.d);
        }
    }
}

void sim_eps() {
    double temp[3] = {15.0, 25.0, 35.0};

    int i, j;
    double eps[4] = {1e-4, 1e-3, 1e-2, 0.1};

    for (j=0; j<3; j++) {
        for (i=0; i<4; i++)
            sim(temp[j], 0.0, 1, eps[i], 0);
    }
}

void sim_sd() {
    double temp[3] = {15.0, 25.0, 35.0};

    int j;

    for (j=0; j<3; j++) {
        sim(temp[j], 4.0, 100, 1e-3, 0);
    }
}

int main(int attr, char *avec[]) {
    spop2_random_init();

    //sim_eps();
    sim_sd();

    return 0;
}
