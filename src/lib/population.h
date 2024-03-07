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
void spop2_random_destroy();

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
#define NOAGE_CONST   9

#define ACC_ARBITER   0
#define AGE_ARBITER   1

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
typedef void (*funqpar)(const number *, const number, double *);

typedef struct arbiter_st *arbiter;
struct arbiter_st {
    parameters fun_pars;
    hazard fun_haz;
    calculator fun_calc;
    stepper fun_step;
    funqpar fun_q_par;
};

arbiter arbiter_init(parameters, hazard, calculator, stepper);

number acc_stepper(number, unsigned int, number);
number age_stepper(number, unsigned int, number);

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

typedef struct population_st *population;
struct population_st {
    size_t nkey;
    char stoch;
    char *types;
    char *arbicodes;
    unsigned int *numpars;
    update fun_update;
    arbiter *arbiters;
    member members;
};

double spop2_version(void);

void spop2_set_eps(double);

void spop2_set_purge(unsigned, double);

typedef void (*transfer)(number *, number, void *);
typedef double (*harvest)(number *);
typedef char (*boolean)(number *);

population spop2_init(char *, char);
void spop2_free(population *);
void spop2_empty(population *);
number spop2_size(population);
number spop2_count(population, boolean);
number spop2_remove(population, number *, double);
number spop2_harvest(population, population, harvest);
char spop2_add(population, number *, number);
char spop2_addpop(population, population);
void spop2_foreach(population, transfer, void *);
void spop2_step(population, double *, number *, number *, population *);
void spop2_print(population);
void spop2_printable(population, int);
unsigned int spop2_buffsize(population);
number *spop2_savestate(population);
population spop2_loadstate_empty(number *);
population spop2_loadstate(number *);

/* ----------------------------------------------------------- *\
 * 
\* ----------------------------------------------------------- */

#endif
