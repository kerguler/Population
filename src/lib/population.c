/*
 * sPop2: a dynamically-structured matrix population model
 * Copyright (C) 2022 Kamil Erguler <k.erguler@cyi.ac.cy>
 * 
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <https://www.gnu.org/licenses/>.
 */

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_cdf.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include "population.h"

/* ----------------------------------------------------------- *\
 * RANDOM
\* ----------------------------------------------------------- */

gsl_rng *RANDOM = 0;
void spop2_random_init() {
    rng_setup("sPop2");
    RANDOM = get_RAND_GSL();
}

/* ----------------------------------------------------------- *\
 * EPS
\* ----------------------------------------------------------- */

double POPSIZE_ROUND_EPS = 1e-13;
void spop2_set_eps(double eps) {
    POPSIZE_ROUND_EPS = eps;
}

double popsize_round(double val) {
    return POPSIZE_ROUND_EPS ? round(val / POPSIZE_ROUND_EPS) * POPSIZE_ROUND_EPS : val;
}

/* ----------------------------------------------------------- *\
\* ----------------------------------------------------------- */

arbiter arbiter_init(parameters fun_pars, hazard fun_haz, calculator fun_calc, stepper fun_step) {
    arbiter step = (arbiter)malloc(sizeof(struct arbiter_st));
    step->fun_pars = fun_pars;
    step->fun_haz = fun_haz;
    step->fun_calc = fun_calc;
    step->fun_step = fun_step;
    return step;
}

void arbiter_free(arbiter *arbiter) {
    free(*arbiter);
}

/* ----------------------------------------------------------- *\
 * Stepper routines
\* ----------------------------------------------------------- */

number acc_stepper(number q, unsigned int d, number k) {
    number q2 = q;
    if (d)
        q2.d = popsize_round(q2.d + (double)d / (double)k.i);
    return q2;
}

number age_stepper(number q, unsigned int d, number k) {
    number q2 = q;
    q2.i++;
    return q2;
}

/* ----------------------------------------------------------- *\
\* ----------------------------------------------------------- */

double acc_hazard_calc(hazard heval, unsigned int d, number q, number k, double theta, const number *qkey) {
    double h0 = d ? heval(d - 1, k, theta) : 0.0;
    double h1 = heval(d, k, theta);
    return h0 == 1.0 ? 1.0 : 1.0 - (1.0 - h1)/(1.0 - h0); // mortality
}

double age_const_calc(hazard heval, unsigned int d, number q, number k, double theta, const number *qkey) {
    return theta; // mortality
}

double age_hazard_calc(hazard heval, unsigned int d, number q, number k, double theta, const number *qkey) {
    double h0 = q.i ? heval(q.i - 1, k, theta) : 0.0;
    double h1 = heval(q.i, k, theta);
    return h0 == 1.0 ? 1.0 : 1.0 - (1.0 - h1)/(1.0 - h0); // mortality
}

double age_custom_calc(hazard heval, unsigned int d, number q, number k, double theta, const number *qkey) {
    return 0.0; // mortality
}

/* ----------------------------------------------------------- *\
\* ----------------------------------------------------------- */

hazpar acc_fixed_pars(double devmn, double devsd) {
    hazpar hz;
    hz.k.i = round(devmn);
    hz.theta = 1.0;
    hz.stay = TRUE;
    return hz;
}

hazpar acc_erlang_pars(double devmn, double devsd) {
    hazpar hz;
    hz.theta = pow(devsd,2) / devmn;
    double kd = devmn / hz.theta;
    if (kd != round(kd)) {
        kd = round(kd);
        hz.theta = devmn / kd;
        // double m = kd * hz.theta;
        // double s = sqrt(hz.theta * m);
        // fprintf(stderr, "Rounding up k to %g to yield mean=%g and sd=%g\n", kd, m, s);
    }
    hz.k.i = (unsigned int)kd;
    hz.stay = TRUE;
    return hz;
}

hazpar acc_pascal_pars(double devmn, double devsd) {
    hazpar hz;
    hz.theta = devmn / pow(devsd,2);
    if (hz.theta > 1.0 || hz.theta < 0.0) {
        fprintf(stderr, "Pascal cannot yield mean=%g and sd=%g\n",devmn,devsd);
        exit(1);
    }
    double kd = devmn * hz.theta / (1.0 - hz.theta);
    if (kd != round(kd)) {
        kd = round(kd);
        hz.theta = kd / (devmn + kd);
    }
    hz.k.i = (unsigned int)kd;
    hz.stay = TRUE;
    return hz;
}

hazpar age_fixed_pars(double devmn, double devsd) {
    hazpar hz;
    hz.k.d = round(devmn);
    hz.theta = 1.0;
    hz.stay = FALSE;
    return hz;
}

hazpar age_const_pars(double devmn, double devsd) {
    hazpar hz;
    hz.k.d = 1.0;
    hz.theta = min(1.0, max(0.0, devmn));
    hz.stay = FALSE;
    return hz;
}

hazpar age_gamma_pars(double devmn, double devsd) {
    hazpar hz;
    hz.theta = pow(devsd,2) / devmn;
    hz.k.d = devmn / hz.theta;
    hz.stay = FALSE;
    return hz;
}

hazpar age_nbinom_pars(double devmn, double devsd) {
    hazpar hz;
    hz.theta = devmn / pow(devsd,2);
    if (hz.theta > 1.0 || hz.theta < 0.0) {
        fprintf(stderr, "Negative binomial cannot yield mean=%g and sd=%g\n",devmn,devsd);
        exit(1);
    }
    hz.k.d = devmn * hz.theta / (1.0 - hz.theta);
    hz.stay = FALSE;
    return hz;
}

hazpar age_custom_pars(double devmn, double devsd) {
    hazpar hz = {.k=numZERO, .theta=1.0, .stay=FALSE};
    return hz;
}

/* ----------------------------------------------------------- *\
\* ----------------------------------------------------------- */

double acc_fixed_haz(unsigned int i, number k, double theta) {
    return (double)((double)i >= theta);
}

double acc_erlang_haz(unsigned int i, number k, double theta) {
    return gsl_cdf_poisson_P(i, ONE/theta);
}

double acc_pascal_haz(unsigned int i, number k, double theta) {
    return 1.0 - pow(theta, i + 1);
}

double age_fixed_haz(unsigned int i, number k, double theta) {
    return (double)((double)i >= k.d);
}

double age_const_haz(unsigned int i, number k, double theta) {
    return theta;
}

double age_gamma_haz(unsigned int i, number k, double theta) {
    // TO DO: The problem of overflow should be solved!
    return gsl_cdf_gamma_P((double)i, k.d, theta);
}

double age_nbinom_haz(unsigned int i, number k, double theta) {
    return i ? gsl_cdf_negative_binomial_P(i-1, theta, k.d) : 0.0;
}

double age_custom_haz(unsigned int i, number k, double theta) {
    return 0.0;
}

/* ----------------------------------------------------------- *\
 * 
\* ----------------------------------------------------------- */

void update_stoch(double p, number *n, number *n2) {
    (*n2).i = gsl_ran_binomial(RANDOM, p, (*n).i);
    (*n).i -= (*n2).i;
}

void update_det(double p, number *n, number *n2) {
    (*n2).d = (*n).d * p;
    (*n).d -= (*n2).d;
}

/* ----------------------------------------------------------- *\
 * sPop2 - Member stack handlers
\* ----------------------------------------------------------- */

member_stack *member_stack_init(unsigned int nkey, char *types, char stoch) {
    member_stack *poptable = (member_stack *)malloc(sizeof(struct member_stack_st));

    poptable->nkey = nkey;

    poptable->keytypes = (char *)malloc(poptable->nkey * sizeof(char));
    poptable->numtype = stoch ? TYP_INT : TYP_FLOAT;

    poptable->key_ids = (unsigned int *)malloc(poptable->nkey * sizeof(unsigned int));

    poptable->key_size = 0;
    poptable->num_size = (stoch ? sizeof(unsigned int) : sizeof(double))/sizeof(char);

    int i;
    for (i=0; i<poptable->nkey; i++) {
        poptable->key_ids[i] = poptable->key_size;
        switch (types[i]) {
            case ACC_ARBITER:
                poptable->key_size += sizeof(double)/sizeof(char);
                poptable->keytypes[i] = TYP_FLOAT;
                break;
            case AGE_ARBITER:
                poptable->key_size += sizeof(unsigned int)/sizeof(char);
                poptable->keytypes[i] = TYP_INT;
                break;
            default:
                fprintf(stderr, "Wrong key type: %d\n", types[i]);
                exit(1);
                break;
        }
    }

    poptable->member_size = poptable->key_size + poptable->num_size;
    poptable->maxmember = MEMBER_BUFF;
    poptable->members = (void *)malloc(poptable->maxmember * poptable->member_size * sizeof(char));
    poptable->nmember = 0;

    return poptable;
}

void member_stack_free(member_stack *poptable) {
    free(poptable->keytypes);
    free(poptable->key_ids);
    free(poptable->members);
    free(poptable);
}

void member_stack_resize(member_stack *poptable) {
    if (poptable->nmember < poptable->maxmember) return;
    if (poptable->maxmember - poptable->nmember < MEMBER_BUFF) return;
    poptable->maxmember = poptable->nmember + MEMBER_BUFF;
    void *tmp = (void *)malloc(poptable->maxmember * poptable->member_size * sizeof(char));
    memcpy(tmp, poptable->members, poptable->nmember * poptable->member_size * sizeof(char));
    free(poptable->members);
    poptable->members = tmp;
}

void *member_stack_search(member_stack *poptable, void *key) {
    void *dst;
    size_t sz = poptable->key_size * sizeof(char);
    int i;
    for (i=0, dst=poptable->members; i<poptable->nmember; i++, dst=(void *)((char *)dst + poptable->member_size))
        if (!memcmp(dst, key, sz))
            return dst;
    return 0;
}

void member_stack_add(member_stack *poptable, number *key_raw, number num) {
    member_stack_resize(poptable);
    //
    void *key = (void *)malloc(poptable->key_size * sizeof(char));
    member_stack_setkey(poptable, key_raw, key);
    void *dst = member_stack_search(poptable, key);
    //
    unsigned int *tmpi;
    double *tmpd;
    //
    if (dst) {
        switch (poptable->numtype) {
            case TYP_INT:
                tmpi = (unsigned int *)((char *)dst + poptable->key_size);
                *tmpi += num.i;
                break;
            case TYP_FLOAT:
                tmpd = (double *)((char *)dst + poptable->key_size);
                *tmpd += num.d;
                break;
            default:
                fprintf(stderr, "Wrong num type: %d\n", poptable->numtype);
                exit(1);
                break;
        }
        free(key);
        return;
    }
    //
    dst = (void *)((char *)(poptable->members) + poptable->nmember * poptable->member_size);
    memcpy(dst, key, poptable->key_size * sizeof(char));
    member_stack_setkey(poptable, key_raw, dst);
    switch (poptable->numtype) {
        case TYP_INT:
            memcpy((char *)dst + poptable->key_size, &num.i, sizeof(unsigned int));
            break;
        case TYP_FLOAT:
            memcpy((char *)dst + poptable->key_size, &num.d, sizeof(double));
            break;
        default:
            fprintf(stderr, "Wrong num type: %d\n", poptable->numtype);
            exit(1);
            break;
    }
    poptable->nmember++;
    //
    free(key);
}

number member_stack_remove(member_stack *poptable, number *key_raw, double frac) {
    void *key = (void *)malloc(poptable->key_size * sizeof(char));
    member_stack_setkey(poptable, key_raw, key);
    void *dst = member_stack_search(poptable, key);
    //
    if (!dst) {
        free(key);
        return numZERO;
    }
    //
    number ret = numZERO;
    void *pnt = (void *)((char *)dst + poptable->key_size);
    //
    unsigned int tmpi;
    double tmpd;
    //
    switch (poptable->numtype) {
        case TYP_INT:
            tmpi = *(unsigned int *)pnt;
            ret.i = gsl_ran_binomial(RANDOM, frac, tmpi);
            tmpi -= ret.i;
            if (tmpi) {
                memcpy((char *)dst + poptable->key_size, &tmpi, sizeof(unsigned int));
            } else if (poptable->nmember > 1) {
                memcpy(dst, (char *)poptable->members + (poptable->nmember - 1) * poptable->member_size, poptable->member_size * sizeof(char));
                poptable->nmember--;
            } else {
                poptable->nmember--;
            }
            break;
        case TYP_FLOAT:
            tmpd = *(double *)pnt;
            ret.d = frac * tmpd;
            tmpd -= ret.d;
            if (tmpd) {
                memcpy((char *)dst + poptable->key_size, &tmpd, sizeof(double));
            } else if (poptable->nmember > 1) {
                memcpy(dst, (char *)poptable->members + (poptable->nmember - 1) * poptable->member_size, poptable->member_size * sizeof(char));
                poptable->nmember--;
            } else {
                poptable->nmember--;
            }
            break;
        default:
            fprintf(stderr, "Wrong num type: %d\n", poptable->numtype);
            exit(1);
            break;
    }
    //
    free(key);
    member_stack_resize(poptable);
    return ret;
}

void member_stack_setkey(member_stack *poptable, number *key_raw, void *dst) {
    int i;
    for (i=0; i<poptable->nkey; i++) {
        switch (poptable->keytypes[i]) {
            case TYP_FLOAT:
                *((double *)((char *)dst + poptable->key_ids[i])) = key_raw[i].d;
                break;
            case TYP_INT:
                *((unsigned int *)((char *)dst + poptable->key_ids[i])) = key_raw[i].i;
                break;
            default:
                fprintf(stderr, "Wrong key type: %d\n", poptable->keytypes[i]);
                exit(1);
                break;
        }
    }
}

void member_stack_getkey_id(member_stack *poptable, void *key, int id, number *key_raw) {
    *key_raw = numZERO;
    switch (poptable->keytypes[id]) {
        case TYP_FLOAT:
            (*key_raw).d = *((double *)((char *)key + poptable->key_ids[id]));
            break;
        case TYP_INT:
            (*key_raw).i = *((unsigned int *)((char *)key + poptable->key_ids[id]));
            break;
        default:
            fprintf(stderr, "Wrong key type: %d\n", poptable->keytypes[id]);
            exit(1);
            break;
    }
}

void member_stack_getkey(member_stack *poptable, void *key, number *key_raw) {
    int i;
    for (i=0; i<poptable->nkey; i++) {
        key_raw[i] = numZERO;
        switch (poptable->keytypes[i]) {
            case TYP_FLOAT:
                key_raw[i].d = *((double *)((char *)key + poptable->key_ids[i]));
                break;
            case TYP_INT:
                key_raw[i].i = *((unsigned int *)((char *)key + poptable->key_ids[i]));
                break;
            default:
                fprintf(stderr, "Wrong key type: %d\n", poptable->keytypes[i]);
                exit(1);
                break;
        }
    }
}

void member_stack_getnum(member_stack *poptable, void *key, number *num) {
    void *tmp = (void *)((char *)key + poptable->key_size);
    *num = numZERO;
    switch (poptable->numtype) {
        case TYP_FLOAT:
            (*num).d = *((double *)tmp);
            break;
        case TYP_INT:
            (*num).i = *((unsigned int *)tmp);
            break;
        default:
            fprintf(stderr, "Wrong key type: %d\n", poptable->numtype);
            exit(1);
            break;
    }
}

void member_stack_numsum(member_stack *poptable, number *numsum) {
    *numsum = numZERO;
    void *dst;
    int i;
    switch (poptable->numtype) {
        case TYP_INT:
            for (i=0, dst=poptable->members; i<poptable->nmember; i++, dst=(void *)((char *)dst + poptable->member_size))
                (*numsum).i += *(unsigned int *)((char *)dst + poptable->key_size);
            break;
        case TYP_FLOAT:
            for (i=0, dst=poptable->members; i<poptable->nmember; i++, dst=(void *)((char *)dst + poptable->member_size))
                (*numsum).d += *(double *)((char *)dst + poptable->key_size);
            break;
        default:
            fprintf(stderr, "Wrong num type: %d\n", poptable->numtype);
            exit(1);
            break;
    }
}

void member_stack_printkey(member_stack *poptable, void *key) {
    void *tmp;
    int i;
    for (i=0, tmp=key; i<poptable->nkey; i++) {
        switch (poptable->keytypes[i]) {
            case TYP_FLOAT:
                printf("%s%g", i ? "," : "", *((double *)tmp));
                tmp = (char *)tmp + sizeof(double)/sizeof(char);
                break;
            case TYP_INT:
                printf("%s%u", i ? "," : "", *((unsigned int *)tmp));
                tmp = (char *)tmp + sizeof(unsigned int)/sizeof(char);
                break;
            default:
                fprintf(stderr, "Wrong key type: %d\n", poptable->keytypes[i]);
                exit(1);
                break;
        }
    }
}

void member_stack_printnum(member_stack *poptable, void *key) {
    void *tmp = (void *)((char *)key + poptable->key_size * sizeof(char));
    switch (poptable->numtype) {
        case TYP_FLOAT:
            printf("%g", *((double *)tmp));
            break;
        case TYP_INT:
            printf("%u", *((unsigned int *)tmp));
            break;
        default:
            fprintf(stderr, "Wrong key type: %d\n", poptable->numtype);
            exit(1);
            break;
    }
}

void member_stack_print(member_stack *poptable) {
    if (!poptable->nmember) {
        printf("Empty population\n");
        return;
    }
    void *dst;
    int i;
    for (i=0, dst = poptable->members; i<poptable->nmember; i++, dst = (void *)((char *)dst + poptable->member_size)) {
        printf("Member{ (");
        member_stack_printkey(poptable, dst);
        switch (poptable->numtype) {
            case TYP_INT:
                printf(") => %u }\n", *(unsigned int *)((char *)dst+poptable->key_size));
                break;
            case TYP_FLOAT:
                printf(") => %g }\n", *(double *)((char *)dst+poptable->key_size));
                break;
            default:
                fprintf(stderr, "Wrong num type: %d\n", poptable->numtype);
                exit(1);
                break;
        }
    }
}

void member_stack_printable(member_stack *poptable, int tm) {
    if (!poptable->nmember) {
        return;
    }
    void *dst;
    int i;
    for (i=0, dst = poptable->members; i<poptable->nmember; i++, dst = (void *)((char *)dst + poptable->member_size)) {
        printf("%d,",tm);
        member_stack_printkey(poptable, dst);
        printf(",");
        member_stack_printnum(poptable, dst);
        printf("\n");
    }
}

/* ----------------------------------------------------------- *\
 * sPop2 - Population handlers
\* ----------------------------------------------------------- */

void spop2_print(population pop) {
    printf("Dynamics: %s\n",pop->stoch?"Stochastic":"Deterministic");
    printf("Key size: %u\n",pop->nkey);
    //
    member_stack_print(pop->poptable);
}

void spop2_printable(population pop, int tm) {
    member_stack_printable(pop->poptable, tm);
}

population spop2_init(char *arbiters, char stoch) {
    int i;
    //
    population pop = (population)malloc(sizeof(struct population_st));
    //
    for (i=0, pop->nkey=0; arbiters[i] != STOP; i++) pop->nkey++;
    if (pop->nkey == 0) {
        fprintf(stderr, "Couldn't find any development time distributions\n");
        exit(1);
    }
    //
    pop->stoch = stoch;
    pop->fun_update = pop->stoch ? update_stoch : update_det;
    //
    pop->types = (char *)calloc(pop->nkey, sizeof(char));
    pop->numpars = (unsigned int *)calloc(pop->nkey, sizeof(unsigned int));
    //
    pop->arbiters = (arbiter *)malloc(pop->nkey * sizeof(struct arbiter_st));
    for (i=0; i < pop->nkey; i++) {
        pop->types[i] = arbiters[i];
        switch (arbiters[i]) {
            case ACC_FIXED:
                pop->arbiters[i] = arbiter_init(acc_fixed_pars, acc_fixed_haz, acc_hazard_calc, acc_stepper);
                pop->types[i] = ACC_ARBITER;
                pop->numpars[i] = 1;
                break;
            case ACC_ERLANG:
                pop->arbiters[i] = arbiter_init(acc_erlang_pars, acc_erlang_haz, acc_hazard_calc, acc_stepper);
                pop->types[i] = ACC_ARBITER;
                pop->numpars[i] = 2;
                break;
            case ACC_PASCAL:
                pop->arbiters[i] = arbiter_init(acc_pascal_pars, acc_pascal_haz, acc_hazard_calc, acc_stepper);
                pop->types[i] = ACC_ARBITER;
                pop->numpars[i] = 2;
                break;
            case AGE_FIXED:
                pop->arbiters[i] = arbiter_init(age_fixed_pars, age_fixed_haz, age_hazard_calc, age_stepper);
                pop->types[i] = AGE_ARBITER;
                pop->numpars[i] = 1;
                break;
            case AGE_CONST:
                pop->arbiters[i] = arbiter_init(age_const_pars, age_const_haz, age_const_calc, age_stepper);
                pop->types[i] = AGE_ARBITER;
                pop->numpars[i] = 1;
                break;
            case AGE_GAMMA:
                pop->arbiters[i] = arbiter_init(age_gamma_pars, age_gamma_haz, age_hazard_calc, age_stepper);
                pop->types[i] = AGE_ARBITER;
                pop->numpars[i] = 2;
                break;
            case AGE_NBINOM:
                pop->arbiters[i] = arbiter_init(age_nbinom_pars, age_nbinom_haz, age_hazard_calc, age_stepper);
                pop->types[i] = AGE_ARBITER;
                pop->numpars[i] = 2;
                break;
            case AGE_CUSTOM:
                pop->arbiters[i] = arbiter_init(age_custom_pars, age_custom_haz, age_custom_calc, age_stepper);
                pop->types[i] = AGE_ARBITER;
                pop->numpars[i] = 0;
                break;
            default:
                fprintf(stderr, "Development time distribution %d not yet implemented\n", arbiters[i]);
                exit(1);
                break;
        }
    }
    //
    pop->poptable = member_stack_init(pop->nkey, pop->types, pop->stoch);
    //
    return pop;
}

void spop2_free(population *pop) {
    int i;
    for (i=0; i < (*pop)->nkey; i++)
        arbiter_free(&((*pop)->arbiters[i]));
    free((*pop)->arbiters);
    free((*pop)->types);
    free((*pop)->numpars);
    member_stack_free((*pop)->poptable);
    spop2_empty(pop);
    free(*pop);
}

void spop2_empty(population *pop) {
    if (!(*pop)) return;
    (*pop)->poptable->nmember = 0;
    member_stack_resize((*pop)->poptable);
}

number spop2_size(population pop) {
    number sz = numZERO;
    member_stack_numsum(pop->poptable, &sz);
    return sz;
}

number spop2_remove(population pop, number *key_raw, double frac) {
    return member_stack_remove(pop->poptable, key_raw, frac);
}

char spop2_add(population pop, number *key_raw, number num) {
    member_stack_add(pop->poptable, key_raw, num);
    return 0;
}

char spop2_addpop(population popdst, population popsrc) {
    if (!popsrc->poptable->nmember) {
        return 0;
    }
    void *src;
    number *key = (number *)malloc(popsrc->nkey * sizeof(number));
    number num;
    int i;
    for (i=0, src = popsrc->poptable->members; i<popsrc->poptable->nmember; i++, src = (void *)((char *)src + popsrc->poptable->member_size)) {
        member_stack_getkey(popsrc->poptable, src, key);
        member_stack_getnum(popsrc->poptable, src, &num);
        if (!memcmp(&num,&numZERO,sizeof(number))) continue;
        member_stack_add(popdst->poptable, key, num);
    }
    //
    free(key);
    return 0;
}

void spop2_foreach(population pop, void (*func)(number *key, number num)) {
    number *key = (number *)malloc(pop->nkey * sizeof(number));
    number num;
    void *dst;
    int i;
    for (i=0, dst = pop->poptable->members; i<pop->poptable->nmember; i++, dst = (void *)((char *)dst + pop->poptable->member_size)) {
        member_stack_getkey(pop->poptable, dst, key);
        member_stack_getnum(pop->poptable, dst, &num);
        func(key, num);        
    }
    free(key);
}

void spop2_step(population pop, double *par, number *survived, number *completed, population *popdone) {
    int i, j;
    //
    hazpar hp;
    unsigned int dev;
    number *q2 = (number *)malloc(pop->nkey * sizeof(number));
    number n2 = numZERO;
    double p;
    number *key = (number *)malloc(pop->nkey * sizeof(number));
    number num;
    member_stack *poptablenext;
    void *dst;
    //
    for (i=0; i < pop->nkey; i++) {
        *survived = numZERO;
        completed[i] = numZERO;
        // 
        switch (pop->numpars[i]) {
            case 0:
                hp = pop->arbiters[i]->fun_pars(0.0, 0.0);
            break;
            case 1:
                hp = pop->arbiters[i]->fun_pars(par[0], 0.0);
                par += 1;
            break;
            case 2:
                hp = pop->arbiters[i]->fun_pars(par[0],par[1]);
                par += 2;
            break;
            default:
                fprintf(stderr, "Hazard functions with more than 2 parameters are not yet supported\n");
                exit(1);
            break;
        }
        // 
        if (!memcmp(&hp,&noHazard,sizeof(struct hazpar_st))) continue;
        //
        poptablenext = member_stack_init(pop->nkey, pop->types, pop->stoch);
        for (j=0, dst = pop->poptable->members; j<pop->poptable->nmember; j++, dst = (void *)((char *)dst + pop->poptable->member_size)) {
            member_stack_getkey(pop->poptable, dst, key);
            member_stack_getnum(pop->poptable, dst, &num);
            if (!memcmp(&num,&numZERO,sizeof(number))) continue;
            //
            if (pop->stoch)
                (*survived).i += num.i;
            else
                (*survived).d += num.d;
            //
            for (dev=0; memcmp(&num,&numZERO,sizeof(number)); ) {
                memcpy(q2, key, pop->nkey * sizeof(number));
                if (pop->arbiters[i]->fun_step)
                    q2[i] = pop->arbiters[i]->fun_step(q2[i], dev, hp.k);
                //
                if (pop->types[i] == ACC_ARBITER ? q2[i].d >= ACCTHR : FALSE) {
                    if (popdone) {
                        spop2_add(popdone[i], q2, num);
                    }
                    if (pop->stoch) {
                        completed[i].i += num.i;
                        (*survived).i -= num.i;
                    } else {
                        completed[i].d += num.d;
                        (*survived).d -= num.d;
                    }
                    num = numZERO;
                } else {
                    p = pop->arbiters[i]->fun_calc(pop->arbiters[i]->fun_haz, dev, q2[i], hp.k, hp.theta, key);
                    pop->fun_update(p, &num, &n2);
                    //
                    if (pop->types[i] == AGE_ARBITER) {
                        if (memcmp(&num, &numZERO, sizeof(number)))
                            member_stack_add(poptablenext, q2, num); // Developing / surviving population
                        if (memcmp(&n2, &numZERO, sizeof(number))) {
                            if (popdone) {
                                spop2_add(popdone[i], q2, n2);
                            }
                            if (pop->stoch) {
                                completed[i].i += n2.i;
                                (*survived).i -= n2.i;
                            } else {
                                completed[i].d += n2.d;
                                (*survived).d -= n2.d;
                            }
                        }
                    } else if (memcmp(&n2,&numZERO,sizeof(number)))
                        member_stack_add(poptablenext, q2, n2); // Developing / surviving population
                    //
                    dev++;
                }
                //
                if (!hp.stay) break;
            }
        }
        //
        if (poptablenext) {
            member_stack_free(pop->poptable);
            pop->poptable = poptablenext;
        }
    }
    //
    free(q2);
    free(key);
}