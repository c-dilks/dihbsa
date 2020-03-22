// test new Modulation class

R__LOAD_LIBRARY(DihBsa)
#include "Modulation.h"
#include "Tools.h"

void testModulationClass() {

  // OPTIONS:
  Int_t pol = Modulation::kUU;
  Bool_t enablePW = false;
  ////////////////

  Modulation * modu;
  Int_t levMax;
  TString titleStr;
  for(int T=2; T<=3; T++) {
    if(pol==Modulation::kUU && T==2) levMax=1;
    else levMax = 0;
    for(int lev=0; lev<=levMax; lev++) {
      titleStr = Form("twist=%d  level=%d",T,lev);
      Tools::PrintTitleBox(titleStr);
      for(int L=0; L<=2; L++) {
        for(int M=-L; M<=L; M++) {
          modu = new Modulation(T,L,M,lev,enablePW,pol);
          printf(" %s\n",(modu->StateTitle()).Data());
          printf("\t%s\n",(modu->GetBaseString()).Data());
          printf("\t%s\n",(modu->Formu()).Data());
        };
        printf("\n");
      };
    };
  };
};
