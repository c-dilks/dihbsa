#include "Binning.h"

ClassImp(Binning)

using namespace std;


Binning::Binning() {
  printf("Instantiating Binning...\n");

  // set up dnp2018 binning
  minIV[vM] = 0;   maxIV[vM] = 3;
  minIV[vX] = 0;   maxIV[vX] = 1.2;
  minIV[vZ] = 0;   maxIV[vZ] = 1.2;
  minIV[vPt] = 0;  maxIV[vPt] = 3;
  for(int v=0; v<nIV; v++) nBins[v]=-1;
  // -- mass
  AddBinBound(vM,minIV[vM]);
  AddBinBound(vM,0.4);
  AddBinBound(vM,0.8);
  AddBinBound(vM,maxIV[vM]);
  // -- x
  AddBinBound(vX,minIV[vX]);
  AddBinBound(vX,0.2);
  AddBinBound(vX,0.4);
  AddBinBound(vX,maxIV[vX]);
  // -- z
  AddBinBound(vZ,minIV[vZ]);
  AddBinBound(vZ,0.4);
  AddBinBound(vZ,0.6);
  AddBinBound(vZ,maxIV[vZ]);
  // -- PhPerp
  AddBinBound(vPt,minIV[vPt]);
  AddBinBound(vPt,0.5);
  AddBinBound(vPt,maxIV[vPt]);

  /*
  if(singleBinMode) {
    for(int v=0; v<nIV; v++) {
      nBins[v]=-1;
      AddBinBound(v,minIV[v]);
      AddBinBound(v,maxIV[v]);
    };
    printf("\n-- SINGLE BIN MODE ENABLED\n");
  };
  */


  IVname[vM] = "M";
  IVname[vX] = "X";
  IVname[vZ] = "Z";
  IVname[vPt] = "Pt";

  IVtitle[vM] = "M_{h}";
  IVtitle[vX] = "x";
  IVtitle[vZ] = "z";
  IVtitle[vPt] = "P_{h}^{perp}";


  PrintBinBounds();


};



void Binning::AddBinBound(Int_t ivIdx, Float_t newBound) {
  if(ivIdx<0 || ivIdx>=nIV) {
    fprintf(stderr,"ERROR: bad Binning::AddBinBound call");
    return;
  };

  bound[ivIdx].push_back(newBound);
  nBins[ivIdx]++;

  return;
};


void Binning::PrintBinBounds() {
  printf("\n");
  for(int v=0; v<nIV; v++) {
    printf("[ %s ] %s bins:  (nbins=%d)\n",IVname[v].Data(),IVtitle[v].Data(),nBins[v]);
    for(int b=0; b<nBins[v]; b++) {
      printf(" bin %d:\t\t%.2f\t%.2f\n",b,bound[v][b],bound[v][b+1]);
    };
  };
  printf("\n");
};


Int_t Binning::GetBin(Int_t ivIdx_, Float_t iv_) {
  if(ivIdx_<0 || ivIdx_>=nIV) {
    fprintf(stderr,"ERROR: bad Binning::GetBin call\n");
    return -1;
  };

  for(int b=0; b<nBins[ivIdx_]; b++) {

    if( iv_ >= bound[ivIdx_].at(b) &&
        iv_ <  bound[ivIdx_].at(b+1) ) {
      return b;
    };
  };

  fprintf(stderr,"ERROR bin not found for %s=%.2f\n",IVname[ivIdx_].Data(),iv_);
  return -1;
};


TString Binning::GetBoundStr(Int_t v_, Int_t b_) {
  TString retStr;
  Float_t lb,ub;
  try {
    lb = bound[v_].at(b_);
    ub = bound[v_].at(b_+1);
  } catch(const std::out_of_range & ex) {
    fprintf(stderr,"ERROR: bad GetBoundStr call\n");
    return "";
  };
  retStr = Form("%s#in[%.2f, %.2f)",IVtitle[v_].Data(),lb,ub);
  return retStr;
};
    

Binning::~Binning() {};

