R__LOAD_LIBRARY(DihBsa)
#include "Asymmetry.h"
#include "Binning.h"
void testStream() {
  TFile * f = new TFile("spinroot/test.root","RECREATE");
  Binning * BS = new Binning(0x34);
  Asymmetry * A = new Asymmetry(BS,0,1,0,0);
  BS->Write("BS");
  A->Write("A");
};

