#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "population.h"

extern gsl_rng *RANDOM;

#define N 10000

/* ----------------------------------------------------------- *\
 * 
\* ----------------------------------------------------------- */

struct key_st {
        int a;
        int b;
};

struct elm {
        struct key_st *key;
        int num;
        UT_hash_handle hmm;
};

void sim1() {
    int i;

    struct key_st *tmp = (struct key_st *)calloc(N, sizeof(struct key_st));
    int *vec = (int *)calloc(N, sizeof(int));

    for (i=0; i<N; i++) {
        struct key_st a;
        a.a = i;
        a.b = i+1;
        tmp[i] = a;
    }

    for (i=0; i<N; i++) 
        vec[i]++;

    free(vec);
}

/* ----------------------------------------------------------- *\
 * 
\* ----------------------------------------------------------- */

void sim2() {
    int i;

    struct elm *tbl = NULL;

    for (i=0; i<N; i++) {
        struct key_st key = {i, i+1};
        struct elm *mem = (struct elm *)malloc(sizeof(struct elm));
        mem->key = &key;
        mem->num = 0;
        HASH_ADD_KEYPTR(hmm, tbl, mem->key, sizeof(struct key_st *), mem);
    }

    struct elm *mm, *tm;
    for (i=0; i<N; i++) {
        HASH_ITER(hmm, tbl, mm, tm) {
            mm->num++;
        }
    }
}

/* ----------------------------------------------------------- *\
 * 
\* ----------------------------------------------------------- */

typedef struct ss_elm_st ss_elm;
struct ss_elm_st {
    number *key;
    number num;
    ss_elm *next;
};

typedef struct ss_list_st ss_list;
struct ss_list_st {
    ss_elm *first;
    ss_elm *last;
};

ss_list ss_list_init() {
    ss_list *list = (ss_list *)malloc(sizeof(ss_list));
    list->first = list->last = 0;
    return *list;
}

char ss_list_add(ss_list *list, number *key, number *num) {
    ss_elm *elm = (ss_elm *)malloc(sizeof(ss_elm));
    elm->key = key;
    elm->num = *num;
    elm->next = 0;
    if (!list->first) {
        list->first = list->last = elm;
    } else {
        list->last->next = elm;
        list->last = elm;
    }
    return 0;
}

void sim3() {
    int i;

    ss_list list = ss_list_init();

    for (i=0; i<N; i++) {
        number key[2] = {{.i=i},{.i=i+1}};
        number num = {.i=0};
        ss_list_add(&list, key, &num);
    }

    ss_elm *elm;
    for (i=0; i<N; i++) {
        for (elm=list.first; elm; elm=elm->next) {
            elm->num.i++;
        }
    }
}

/* ----------------------------------------------------------- *\
 * 
\* ----------------------------------------------------------- */

typedef struct aa_elm_st aa_elm;
struct aa_elm_st {
    number *key;
    number num;
};

aa_elm aa_elm_init(number *key, number *num) {
    aa_elm elm;
    elm.key = key;
    elm.num = *num;
    return elm;
}

void sim4() {
    int i;

    aa_elm *list = (aa_elm *)malloc(N * sizeof(struct aa_elm_st));
    
    for (i=0; i<N; i++) {
        number key[2] = {{.i=i},{.i=i+1}};
        number num = {.i=0};
        list[i] = aa_elm_init(key, &num);
    }

    for (i=0; i<N; i++) {
        list[i].num.i++;
    }

    for (i=0; i<N; i++) {
        list[i].key[1].i++;
    }
}

/* ----------------------------------------------------------- *\
 * 
\* ----------------------------------------------------------- */

void sim5() {
    char arbiters[4] = {AGE_CONST, STOP};
    population pop = spop2_init(arbiters, DETERMINISTIC);

    int i;
    for (i=0; i<N; i++) {
        number key_raw[1] = {{.i=i}};
        number num = {.d=0.0};
        member_stack_add(pop->poptable, key_raw, num);
    }

    for (i=0; i<N; i++) {
        double *dst = (double *)((char *)pop->poptable->members + i * pop->poptable->member_size + pop->poptable->key_size);
        (*dst)++;
    }
}

/* ----------------------------------------------------------- *\
 * 
\* ----------------------------------------------------------- */

void sim6() {
    char arbiters[4] = {AGE_CONST, ACC_ERLANG, AGE_CUSTOM, STOP};
    population pop = spop2_init(arbiters, DETERMINISTIC);

    {
    number key_raw[3] = {{.i=10},{.d=0.4},{.i=3}};
    number num = {.d=0.7};
    member_stack_add(pop->poptable, key_raw, num);
    }

    {
    number key_raw[3] = {{.i=11},{.d=0.5},{.i=6}};
    number num = {.d=1.7};
    member_stack_add(pop->poptable, key_raw, num);
    }

    {
    number key_raw[3] = {{.i=10},{.d=0.4},{.i=3}};
    number num = {.d=0.9};
    member_stack_add(pop->poptable, key_raw, num);
    }

    {
    number key_raw[3] = {{.i=10},{.d=0.4},{.i=3}};
    number num = member_stack_remove(pop->poptable, key_raw, 1.0);
    printf("Removed: %g\n", num.d);
    }

    {
    number key_raw[3] = {{.i=11},{.d=0.4},{.i=6}};
    number num = member_stack_remove(pop->poptable, key_raw, 1.0);
    printf("Removed: %g\n", num.d);
    }

    member_stack_print(pop->poptable);
}

/* ----------------------------------------------------------- *\
 * 
\* ----------------------------------------------------------- */

int main(int attr, char *avec[]) {
    spop2_random_init();

    sim6();

    return 0;
}
