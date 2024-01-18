#include <stdio.h>
#include <stdlib.h>
#include "population.h"

void run_full() {
    char arbiters[3] = {ACC_ERLANG, ACC_ERLANG, STOP};
    population pop1a = spop2_init(arbiters, DETERMINISTIC);
    population pop1b = spop2_init(arbiters, DETERMINISTIC);

    number key[2] = {numZERO,numZERO};
    number num = {.d=10.0};
    spop2_add(pop1a, key, num);
    spop2_add(pop1b, key, num);

    printf("0,%d,%g,%g,%g,%g\n",0,spop2_size(pop1a).d,0.0,spop2_size(pop1b).d,0.0);

    double par[8] = {20.0, 10.0, 10.0, 5.0,
                     10.0,  5.0,  5.0, 2.5};

    number sz1a, cm1a, sz1b, cm1b;
    int i;
    for (i=1; i<20; i++) {
        spop2_step(pop1a, par, &sz1a, &cm1a, 0);
        spop2_step(pop1b, par+4, &sz1b, &cm1b, 0);
        printf("0,%d,%g,%g,%g,%g\n",i,sz1a.d,cm1a.d,sz1b.d,cm1b.d);
    }
}

void run_part() {
    char arbiters[3] = {ACC_ERLANG, ACC_ERLANG, STOP};
    population pop = spop2_init(arbiters, DETERMINISTIC);

    number key[2] = {numZERO,numZERO};
    number num = {.d=10.0};
    spop2_add(pop, key, num);

    printf("%d,%g,%g\n",0,spop2_size(pop).d,0.0);

    double par[4] = {20.0, 10.0, 10.0, 5.0};

    number sz, cm;
    int i;
    for (i=1; i<20; i++) {
        spop2_step(pop, par, &sz, &cm, 0);
        printf("1,%d,%g,%g\n",i,sz.d,cm.d);
    }

    printf("\n");
    spop2_printable(pop,0);

    number *state = spop2_savestate(pop);
    printf("\nBuffer size: %u bytes\n", spop2_buffsize(pop));

    population pop2 = spop2_loadstate(state);

    printf("\n");
    spop2_printable(pop2,0);
}

void run_save() {
    char arbiters[3] = {ACC_ERLANG, ACC_ERLANG, STOP};
    population pop1a = spop2_init(arbiters, DETERMINISTIC);
    population pop1b = spop2_init(arbiters, DETERMINISTIC);

    number key[2] = {numZERO,numZERO};
    number num = {.d=10.0};
    spop2_add(pop1a, key, num);
    spop2_add(pop1b, key, num);

    printf("1,%d,%g,%g,%g,%g\n",0,spop2_size(pop1a).d,0.0,spop2_size(pop1b).d,0.0);

    double par[8] = {20.0, 10.0, 10.0, 5.0,
                     10.0,  5.0,  5.0, 2.5};

    number sz1a, cm1a, sz1b, cm1b;
    int i;
    for (i=1; i<10; i++) {
        spop2_step(pop1a, par, &sz1a, &cm1a, 0);
        spop2_step(pop1b, par+4, &sz1b, &cm1b, 0);
        printf("1,%d,%g,%g,%g,%g\n",i,sz1a.d,cm1a.d,sz1b.d,cm1b.d);
    }

    //
    FILE *file;
    number *buff = 0;
    unsigned int buffsz;
    //
    file = fopen("test.bin","wb");
    rewind(file);
    //
    buffsz = spop2_buffsize(pop1a);
    buff = spop2_savestate(pop1a);
    fwrite(&buffsz, sizeof(unsigned int), 1, file);
    fwrite(buff, buffsz, 1, file);
    free(buff);
    //
    buffsz = spop2_buffsize(pop1b);
    buff = spop2_savestate(pop1b);
    fwrite(&buffsz, sizeof(unsigned int), 1, file);
    fwrite(buff, buffsz, 1, file);
    free(buff);
    //
    fclose(file);
    //
    file = fopen("test.bin","rb");
    rewind(file);
    //
    fread(&buffsz, sizeof(unsigned int), 1, file);
    buff = (number *)malloc(buffsz);
    fread(buff, buffsz, 1, file);
    population pop2a = spop2_loadstate(buff);
    free(buff);
    //
    fread(&buffsz, sizeof(unsigned int), 1, file);
    buff = (number *)malloc(buffsz);
    fread(buff, buffsz, 1, file);
    population pop2b = spop2_loadstate(buff);
    free(buff);
    //
    fclose(file);
    //

    for (i=10; i<20; i++) {
        spop2_step(pop2a, par, &sz1a, &cm1a, 0);
        spop2_step(pop2b, par+4, &sz1b, &cm1b, 0);
        printf("1,%d,%g,%g,%g,%g\n",i,sz1a.d,cm1a.d,sz1b.d,cm1b.d);
    }
}

int main(void) {
    run_full();
    run_save();

    return 0;
}