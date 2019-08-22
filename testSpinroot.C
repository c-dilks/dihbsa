R__LOAD_LIBRARY(DihBsa)
#include "Asymmetry.h"
void testSpinroot() {
  TFile * infile = new TFile("spinroot/spin.4037_4.root","READ");
  infile->ls();
  Asymmetry * a = (Asymmetry*) infile->Get("A_M0");
  a->PrintSettings();
};

