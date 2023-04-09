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

number numACCTHR = {.d=1.0};

/* ----------------------------------------------------------- *\
 * RANDOM
\* ----------------------------------------------------------- */

gsl_rng *RANDOM = 0;
void spop2_random_init() {
    rng_setup("sPop2");
    RANDOM = get_RAND_GSL();
}

void spop2_random_destroy() {
    rng_destroy();
    RANDOM = 0;
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
 * Member key handlers
\* ----------------------------------------------------------- */

number *key_init(number *key_raw, unsigned int nkey, char *types) {
    number *key = (number *)calloc(nkey, sizeof(number));
    int i = 0;
    for (; i < nkey; i++) {
        key[i] = key_raw[i];
        if (types[i] == ACC_ARBITER)
            key[i].d = popsize_round(key_raw[i].d);
    }
    return key;
}

void key_add(member *tbl, number *key, number num, unsigned int nkey, char stoch) {
    member qnt = NULL;
    size_t sz = nkey * sizeof(number);
    if (tbl && *tbl != NULL) {
        HASH_FIND(hh, *tbl, key, sz, qnt);
    }
    if (tbl && *tbl != NULL && qnt) {
        if (stoch)
            qnt->num.i += num.i;
        else
            qnt->num.d += num.d;
    } else {
        member elm = (member)calloc(1, sizeof(struct member_st));
        elm->key = (number *)calloc(nkey, sizeof(number));
        memcpy(elm->key, key, sz);
        elm->num = num;
        HASH_ADD_KEYPTR(hh, *tbl, elm->key, sz, elm);
    }
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
 * sPop2 - Population handlers
\* ----------------------------------------------------------- */

void spop2_print(population pop) {
    printf("Dynamics: %s\n",pop->stoch?"Stochastic":"Deterministic");
    printf("Key size: %zu\n",pop->nkey);
    if (!(pop->members)) {
        printf("Empty population\n");
        return;
    }
    member elm, tmp;
    HASH_ITER(hh, pop->members, elm, tmp) {
        printf("Member{ (");
        int i;
        for (i = 0; i < pop->nkey; i++)
            if (pop->types[i] == ACC_ARBITER)
                printf("%s%g", i ? "," : "", elm->key[i].d);
            else
                printf("%s%u", i ? "," : "", elm->key[i].i);
        if (pop->stoch)
            printf(") => %u }\n",elm->num.i);
        else
            printf(") => %g }\n",elm->num.d);
    }
}


void spop2_printable(population pop, int tm) {
    member elm, tmp;
    HASH_ITER(hh, pop->members, elm, tmp) {
        printf("%d,", tm);
        int i;
        for (i = 0; i < pop->nkey; i++)
            if (pop->types[i] == ACC_ARBITER)
                printf("%s%g", i ? "," : "", elm->key[i].d);
            else
                printf("%s%u", i ? "," : "", elm->key[i].i);
        if (pop->stoch)
            printf(",%u\n",elm->num.i);
        else
            printf(",%g\n",elm->num.d);
    }
}

population spop2_init(char *arbiters, char stoch) {
    int i;
    //
    population pop = (population)malloc(sizeof(struct population_st));
    //
    for (i=0, pop->nkey=0; arbiters[i] != STOP; i++) pop->nkey++;
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
            case NOAGE_CONST:
                pop->arbiters[i] = arbiter_init(age_const_pars, age_const_haz, age_const_calc, 0);
                pop->types[i] = AGE_ARBITER;
                pop->numpars[i] = 1;
                break;
            default:
                fprintf(stderr, "Development time distribution %d not yet implemented\n", arbiters[i]);
                exit(1);
                break;
        }
    }
    //
    pop->members = NULL;
    return pop;
}

void spop2_free(population *pop) {
    int i;
    for (i=0; i < (*pop)->nkey; i++)
        arbiter_free(&((*pop)->arbiters[i]));
    free((*pop)->arbiters);
    free((*pop)->types);
    free((*pop)->numpars);
    spop2_empty(pop);
    free(*pop);
}

void spop2_empty(population *pop) {
    if (!(*pop)) return;
    member elm, tmp;
    HASH_ITER(hh, (*pop)->members, elm, tmp) {
        HASH_DEL((*pop)->members, elm);
        member_free(elm);
    }
}

number spop2_size(population pop) {
    number sz = numZERO;
    member elm, tmp;
    if (pop->stoch) {
        HASH_ITER(hh, pop->members, elm, tmp) {
            sz.i += elm->num.i;
        }    
    } else {
        HASH_ITER(hh, pop->members, elm, tmp) {
            sz.d += elm->num.d;
        }    
    }
    return sz;
}

number spop2_remove(population pop, number *key, double frac) {
    member qnt = NULL;
    HASH_FIND(hh, pop->members, key, pop->nkey * sizeof(number), qnt);
    if (qnt != NULL) {
        number ret = qnt->num;
        //
        if (pop->stoch) {
            ret.i = gsl_ran_binomial(RANDOM, frac, ret.i);
            qnt->num.i -= ret.i;
        } else {
            ret.d *= frac;
            qnt->num.d -= ret.d;
        }
        //
        if (!memcmp(&(qnt->num), &numZERO, sizeof(number))) {
            HASH_DEL(pop->members, qnt);
            member_free(qnt);
        }
        //
        return ret;
    }
    return numZERO;
}

char spop2_add(population pop, number *key_raw, number num) {
    if (!((pop->stoch && num.i > 0) || (!pop->stoch && num.d > 0.0))) return 0;
    //
    number *key = key_init(key_raw, pop->nkey, pop->types);
    key_add(&(pop->members), key, num, pop->nkey, pop->stoch);
    //
    return 0;
}

char spop2_addpop(population popto, population popfrom) {
    member elm, tmp;
    HASH_ITER(hh, popfrom->members, elm, tmp) {
        spop2_add(popto, elm->key, elm->num);
    }    
    //
    return 0;
}

void spop2_foreach(population pop, transfer fun, void *opt) {
    member elm, tmp;
    HASH_ITER(hh, pop->members, elm, tmp) {
        fun(elm->key, elm->num, opt);
    }    
}

void spop2_step(population pop, double *par, number *survived, number *completed, population *popdone) {
    int i;
    //
    size_t sz = pop->nkey * sizeof(number);
    hazpar hp;
    member elm = NULL, tmp = NULL, poptablenext = NULL;
    number *q2 = (number *)malloc(sz);
    number n2 = numZERO;
    unsigned int dev;
    double p;
    //
    for (i=0; i<pop->nkey; i++) completed[i] = numZERO;
    //
    for (i=0; i < pop->nkey; i++) {
        (*survived) = numZERO;
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
        poptablenext = NULL;
        HASH_ITER(hh, pop->members, elm, tmp) {
            if (!memcmp(&(elm->num),&numZERO,sizeof(number))) continue;
            if (pop->stoch)
                (*survived).i += elm->num.i;
            else
                (*survived).d += elm->num.d;
            //
            for (dev=0; memcmp(&(elm->num),&numZERO,sizeof(number)); ) {
                memcpy(q2, elm->key, sz);
                if (pop->arbiters[i]->fun_step)
                    q2[i] = (pop->types[i] == ACC_ARBITER) && (hp.k.i < 1) ? numACCTHR : pop->arbiters[i]->fun_step(q2[i], dev, hp.k);
                //
                if (pop->types[i] == ACC_ARBITER ? q2[i].d >= ACCTHR : FALSE) {
                    if (popdone) {
                        spop2_add(popdone[i], q2, elm->num);
                    }
                    if (pop->stoch) {
                        completed[i].i += elm->num.i;
                        (*survived).i -= elm->num.i;
                    } else {
                        completed[i].d += elm->num.d;
                        (*survived).d -= elm->num.d;
                    }
                    elm->num = numZERO;
                } else {
                    p = pop->arbiters[i]->fun_calc(pop->arbiters[i]->fun_haz, dev, q2[i], hp.k, hp.theta, elm->key);
                    pop->fun_update(p, &(elm->num), &n2);
                    //
                    if (pop->types[i] == AGE_ARBITER) {
                        if (memcmp(&(elm->num), &numZERO, sizeof(number)))
                            key_add(&poptablenext, q2, elm->num, pop->nkey, pop->stoch); // Developing / surviving population
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
                        key_add(&poptablenext, q2, n2, pop->nkey, pop->stoch); // Developing / surviving population
                    //
                    dev++;
                }
                //
                if (!hp.stay) break;
            }
        }
        //
        if (poptablenext) {
            if (pop->members) spop2_empty(&pop);
            pop->members = poptablenext;
        }
    }
    //
    free(q2);
}