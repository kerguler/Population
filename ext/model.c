#include "population.h"

#define NumPar 3
#define NumMet 1
#define NumEnv 1

void init(int *ne, int *np, int *nm) {
  *ne = NumEnv;
  *np = NumPar;
  *nm = NumMet;
}

void destroy(void) {
    
}

void parameters(char **names) {
  char temp[NumMet+NumPar][256] = {
    "species",
    "init","mean","std"
  };
  int i;
  for (i=0; i<(NumMet+NumPar); i++)
    names[i] = strdup(temp[i]);
}

void sim(double               *envar,
         double               *param,
         int                 *finalT,
         double              *result,
         int                *success) {
  double *temperature = envar + 0*(*finalT);
  
}
