R__LOAD_LIBRARY(DihBsa)
#include "Asymmetry.h"
void testSpinroot() {
  TFile * infile = new TFile("spinroot/spin.4039_4.root","READ");
  /*
  infile->ls();
  Asymmetry * a = (Asymmetry*) infile->Get("A_M0");
  a->PrintSettings();
  */
  TIter i(gDirectory->GetListOfKeys());
  TKey * k;
  Asymmetry * a;
  while(( k = (TKey*)i() )) {
    printf("name=%s  className=%s\n", k->GetName(), k->GetClassName() );
    printf(" sizeof = %d\n",k->Sizeof());
    printf(" nbytes = %d\n",k->GetNbytes());
    //printf(" @ %p\n",(void*)k->ReadObj());
    a = (Asymmetry*) infile->Get("A_M0");
    a->SetDirectory(gROOT);
  };
};

