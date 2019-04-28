#include "Asymmetry.h"

ClassImp(Asymmetry)

using namespace std;


Asymmetry::Asymmetry() {
  printf("Instantiating Asymmetry...\n");

  // set up (default) binning
  for(int v=0; v<nIV; v++) nBins[v]=-1;
  // -- mass
  AddBinBound(vM,0);
  AddBinBound(vM,0.5);
  AddBinBound(vM,1);
  AddBinBound(vM,3);
  // -- x
  AddBinBound(vX,0);
  AddBinBound(vX,0.5);
  AddBinBound(vX,1);
  // -- z
  AddBinBound(vZ,0);
  AddBinBound(vZ,0.5);
  AddBinBound(vZ,1);


  IVname[vM] = "M";
  IVname[vX] = "X";
  IVname[vZ] = "Z";

  IVtitle[vM] = "M_{h}";
  IVtitle[vX] = "x";
  IVtitle[vZ] = "z";

  PrintBinBounds();


  // instantiate 1-d distributions
  //for(int v=0; v<nIV; v++) {
};



void Asymmetry::AddBinBound(Int_t iv, Float_t newBound) {
  if(iv<0 || iv>=nIV) {
    fprintf(stderr,"ERROR: bad Asymmetry::AddBinBound call");
    return;
  };

  bound[iv][++nBins[iv]] = newBound;

  return;
};


void Asymmetry::PrintBinBounds() {
  for(int v=0; v<nIV; v++) {
    printf("[---] %s bins:  (nbins=%d)\n",IVtitle[v].Data(),nBins[v]);
    for(int b=0; b<nBins[v]; b++) {
      printf(" bin %d:\t\t%.2f\t%.2f\n",b,bound[v][b],bound[v][b+1]);
    };
  };
};

    



Asymmetry::~Asymmetry() {
};

