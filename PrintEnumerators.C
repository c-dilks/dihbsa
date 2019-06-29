R__LOAD_LIBRARY(src/DihBsa)

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
  Int_t pa,pb;
  TString outf = "pairs.list";
  TString outft = outf+".tmp";
  gSystem->RedirectOutput(outft,"w");
  for(int o1=0; o1<nObservables; o1++) {
    for(int o2=0; o2<nObservables; o2++) {
      if(IterPair(o1,o2,pa,pb)) {
        printf("0x%x %s\n", EncodePairType(pa,pb), PairName(pa,pb).Data() );
      };
    };
  };
  gSystem->RedirectOutput(0);
  gROOT->ProcessLine(TString(".! sort " + outft + " | uniq > " + outf));
  gROOT->ProcessLine(TString(".! rm " + outft));
  gROOT->ProcessLine(TString(".! cat " + outf));
};

