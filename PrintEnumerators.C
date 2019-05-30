R__LOAD_LIBRARY(DihBsa)

#include "Constants.h"

void PrintEnumerators() {
  printf("\n");

  // print particles
  printf("particles\n");
  printf("---------\n");
  for(int p=0; p<nParticles; p++) {
    printf("%d = %s  (pid=%d)\n",p,PartName(p).Data(),PartPID(p));
  };
  printf("\n");

  // print pair types
  printf("dihadron pair types\n");
  printf("-------------------\n");
  for(int p=0; p<nPairType; p++) {
    printf("%d = %s  (%s)\n",p,pairTitle(p).Data(),pairName(p).Data());
  };
  printf("\n");


};

