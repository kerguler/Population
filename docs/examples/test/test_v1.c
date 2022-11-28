#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "population.h"

extern gsl_rng *RANDOM;

#define N 10000

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

int main(int attr, char *avec[]) {
    spop2_random_init();

    sim4();

    return 0;
}