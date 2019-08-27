R__LOAD_LIBRARY(DihBsa)
#include "Constants.h"
#include "Tools.h"

void Test() {

  // test PairSame
  /*
  Int_t m = 2;
  for(int i=0; i<m; i++) {
    for(int j=0; j<m; j++) {
      for(int k=0; k<m; k++) {
        for(int l=0; l<m; l++) {
          printf("%d %d =?= %d %d -- %s\n",i,j,k,l,Tools::PairSame(i,j,k,l)?"yes":"no");
        };
      };
    };
  };
  */

  // test dihHadIdx
  for(int i=0; i<nParticles; i++) {
    for(int j=0; j<nParticles; j++) {
      printf("%s %s --> %s\n",
        PartName(i).Data(),
        PartName(j).Data(),
        PairName(i,j).Data()
      );
      /*
      printf("%s %s --> %s\n",
        PartTitle(i).Data(),
        PartTitle(j).Data(),
        PairTitle(i,j).Data()
      );
      */
    };
  };
};
