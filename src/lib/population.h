#ifndef POPULATION_H
#define POPULATION_H

/* ----------------------------------------------------------- *\
 * 
\* ----------------------------------------------------------- */

#define max(a,b) ((a)>(b)?(a):(b))
#define min(a,b) ((a)<(b)?(a):(b))

/* ----------------------------------------------------------- *\
 * 
\* ----------------------------------------------------------- */

#include "uthash.h"

/* ----------------------------------------------------------- *\
 * 
\* ----------------------------------------------------------- */

#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
gsl_rng *get_RAND_GSL(void);
void rng_setup(char*);
void rng_setup_seed(unsigned int, char*);
void rng_destroy(void);
double rng_exponential(const double);
void spop2_random_init();

/* ----------------------------------------------------------- *\
 * 
\* ----------------------------------------------------------- */

#define EPS_GAMMA 1e-14

#define n_MAX 400.0
#define n_STEP 1.0
#define n_i_MAX 400
#define row_SIZE 400
#define mean_MAX 201.0
#define mean_STEP 0.01
#define mean_i_MAX 20100
#define gamma_SIZE 8040000

#define MAX_DAYS 1000

#define gamma_matrix_sd 0.375

void set_gamma_mem(uint64_t);
void gamma_dist_destroy(void);
void gamma_dist_check(void);
double gamma_pdf(double, double, double);
char gamma_dist_hash(double, double, double, double *);
double gamma_dist_prob(double, double, double);
double gamma_dist_matrix(double, double);
void prepare_gamma_matrix(void);
double nbinom_prob(unsigned int, double, double);
double nbinom_dist_prob(double, double, unsigned int);

/* ----------------------------------------------------------- *\
 * 
\* ----------------------------------------------------------- */

#define TRUE          1
#define FALSE         0

#define ONE           1.0
#define ZERO          0.0

#define STOCHASTIC    1
#define DETERMINISTIC 0

#define ACCTHR        1.0

#define STOP          0
#define ACC_FIXED     1
#define ACC_ERLANG    2
#define ACC_PASCAL    3
#define AGE_FIXED     4
#define AGE_CONST     5
#define AGE_GAMMA     6
#define AGE_NBINOM    7
#define AGE_CUSTOM    8

#define ACC_ARBITER   0
#define AGE_ARBITER   1

#define TYP_INT       0
#define TYP_FLOAT     1

/* ----------------------------------------------------------- *\
 * 
\* ----------------------------------------------------------- */

/*
 * Data type for abundance:
 *  unsigned int: stochastic simulations
 *  double:       deterministic simulations
 */
typedef union {
  unsigned int i;
  double d;
} number;

static number numZERO;

/* ----------------------------------------------------------- *\
 * 
\* ----------------------------------------------------------- */

typedef struct hazpar_st hazpar;
struct hazpar_st {
    number k;
    double theta;
    char stay;
};

static hazpar noHazard;

/* ----------------------------------------------------------- *\
 * 
\* ----------------------------------------------------------- */

typedef hazpar (*parameters)(double, double);
typedef double (*hazard)(unsigned int, number, double);
typedef double (*calculator)(hazard, unsigned int, number, number, double, const number *);
typedef number (*stepper)(number, unsigned int, number);

typedef struct arbiter_st *arbiter;
struct arbiter_st {
    parameters fun_pars;
    hazard fun_haz;
    calculator fun_calc;
    stepper fun_step;
};

arbiter arbiter_init(parameters, hazard, calculator, stepper);

number acc_stepper(number, unsigned int, number);

double acc_hazard_calc(hazard, unsigned int, number, number, double, const number *);
double age_const_calc(hazard, unsigned int, number, number, double, const number *);
double age_hazard_calc(hazard, unsigned int, number, number, double, const number *);
double age_custom_calc(hazard, unsigned int, number, number, double, const number *);

hazpar acc_fixed_pars(double, double);
hazpar acc_erlang_pars(double, double);
hazpar acc_pascal_pars(double, double);
hazpar age_fixed_pars(double, double);
hazpar age_const_pars(double, double);
hazpar age_gamma_pars(double, double);
hazpar age_nbinom_pars(double, double);
hazpar age_custom_pars(double, double);

double acc_fixed_haz(unsigned int, number, double);
double acc_erlang_haz(unsigned int, number, double);
double acc_pascal_haz(unsigned int, number, double);
double age_fixed_haz(unsigned int, number, double);
double age_const_haz(unsigned int, number, double);
double age_gamma_haz(unsigned int, number, double);
double age_nbinom_haz(unsigned int, number, double);
double age_custom_haz(unsigned int, number, double);

/* ----------------------------------------------------------- *\
 * 
\* ----------------------------------------------------------- */

typedef struct member_st *member;
struct member_st {
    number *key;
    number num;
    UT_hash_handle hh;
};

/* ----------------------------------------------------------- *\
 * 
\* ----------------------------------------------------------- */

#define key_free(x) {free(x);}
#define member_free(x) {key_free(x->key); free(x);}
number *key_init(number *, unsigned int, char *);
void key_add(member *, number *, number, unsigned int, char);

/* ----------------------------------------------------------- *\
 * 
\* ----------------------------------------------------------- */

typedef void (*update)(double, number *, number *);
void update_stoch(double p, number *, number *);
void update_det(double p, number *, number *);

/* ----------------------------------------------------------- *\
 * 
\* ----------------------------------------------------------- */

#define MEMBER_BUFF 100

typedef struct member_stack_st member_stack;
struct member_stack_st {
    unsigned int nkey;
    char *keytypes;
    char numtype;
    unsigned int key_size;
    unsigned int num_size;
    unsigned int member_size;
    unsigned int maxmember;
    unsigned int nmember;
    void *members;
};

member_stack *member_stack_init(unsigned int, char *, char);
void member_stack_free(member_stack *);
void member_stack_resize(member_stack *);
void member_stack_setkey(member_stack *, number *, void *);
void member_stack_getkey(member_stack *, void *, number *);
void *member_stack_search(member_stack *, void *);
void member_stack_add(member_stack *, number *, number);
number member_stack_remove(member_stack *, number *, double);
void member_stack_numsum(member_stack *, number *);
void member_stack_printkey(member_stack *, void *);
void member_stack_printnum(member_stack *, void *);
void member_stack_print(member_stack *);
void member_stack_printable(member_stack *);

/* ----------------------------------------------------------- *\
 * 
\* ----------------------------------------------------------- */

typedef struct population_st *population;
struct population_st {
    unsigned int nkey;
    char stoch;
    char *types;
    unsigned int *numpars;
    update fun_update;
    arbiter *arbiters;
    member_stack *poptable;
    member members;
};

void spop2_set_eps(double);

population spop2_init(char *, char);
void spop2_free(population *);
void spop2_empty(population *);
number spop2_size(population);
number spop2_remove(population, number *, double);
char spop2_add(population, number *, number);
void spop2_step(population, double *, number *, number *, population *);
void spop2_print(population);
void spop2_printable(population, int);

/* ----------------------------------------------------------- *\
 * 
\* ----------------------------------------------------------- */

#endif
