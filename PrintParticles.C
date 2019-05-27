R__LOAD_LIBRARY(DihBsa)

#include "Constants.h"

void PrintParticles() {
  for(int p=0; p<nParticles; p++) {
    printf("%d = %s  (pid=%d)\n",p,PartName(p).Data(),PartPID(p));
  };
};

