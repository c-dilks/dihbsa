// test new Modulation class

R__LOAD_LIBRARY(DihBsa)
#include "Modulation.h"

void testModulationClass() {
  Modulation * m = new Modulation();
  //m->enableTheta = false;
  for(int T=2; T<=3; T++) {
    printf("\n twist %d:\n",T);
    for(int L=0; L<=2; L++) {
      for(int M=-L; M<=L; M++) {
        printf(" |%d,%d>: \t%s\n",L,M,(m->BaseString(T,L,M)).Data());
        printf("          \t%s\n",(m->BuildTF3formu(T,L,M)).Data());
      };
      printf("\n");
    };
  };
};
