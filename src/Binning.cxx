#include "Binning.h"

ClassImp(Binning)

using namespace std;


Binning::Binning(Int_t pairType_) {
  //printf("Instantiating Binning...\n");

  // get hadron indices (bin bounds depends on hadron type)
  DecodePairType(pairType_,whichHad[qA],whichHad[qB]);
  numKaons = 0;
  for(int h=0; h<2; h++) { if(whichHad[h]==kKp || whichHad[h]==kKm) numKaons++; };

  // set up dnp2018 binning
  minIV[vM] = 0;   maxIV[vM] = 3;
  minIV[vX] = 0;   maxIV[vX] = 1;
  minIV[vZ] = 0;   maxIV[vZ] = 1;
  minIV[vPt] = 0;  maxIV[vPt] = 3;
  for(int v=0; v<nIV; v++) nBins[v]=-1;
  // -- mass
  AddBinBound(vM,minIV[vM]);
  if(numKaons==0) {
    /*
    AddBinBound(vM,0.4);
    AddBinBound(vM,0.8);
    */
    /*
    AddBinBound(vM,0.4);
    AddBinBound(vM,0.6);
    AddBinBound(vM,0.8);
    AddBinBound(vM,1.0);
    AddBinBound(vM,1.2);
    AddBinBound(vM,1.5);
    */
    ///*
    AddBinBound(vM,0.77); // rho mass
    //*/
  } else if(numKaons==1) {
    AddBinBound(vM,0.85);
    AddBinBound(vM,1.1);
  } else {
    AddBinBound(vM,1.25);
  };
  AddBinBound(vM,maxIV[vM]);
  // -- x
  AddBinBound(vX,minIV[vX]);
  if(numKaons==0) {
    AddBinBound(vX,0.2);
    AddBinBound(vX,0.4);
  } else {
    AddBinBound(vX,0.2);
    AddBinBound(vX,0.3);
  };
  AddBinBound(vX,maxIV[vX]);
  // -- z
  AddBinBound(vZ,minIV[vZ]);
  if(numKaons==0) {
    AddBinBound(vZ,0.4);
    AddBinBound(vZ,0.6);
  } else {
    AddBinBound(vZ,0.5);
    AddBinBound(vZ,0.7);
  };
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


  //PrintBinBounds();


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
      printf(" bin %d:\t\t%.2f\t%.2f\n",b,bound[v].at(b),bound[v].at(b+1));
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


Int_t Binning::GetColor(Int_t v_) {
  switch(v_) {
    //case vM: return kRed;
    //case vX: return kGreen+1;
    //case vZ: return kViolet+2;
    //case vPt: return kAzure;
    default: return kBlack;
  };
};


Float_t Binning::GetAziMax(Int_t v_, Int_t b_) {

  if(v_==vZ) {
    if(b_==0) return 3.950000;
    if(b_==1) return 4.450000;
    if(b_==2) return 4.350000;
  };

  if(v_==vX) {
    if(b_==0) return 4.350000;
    if(b_==1) return 4.450000;
    if(b_==2) return 3.450000;
  };

  if(v_==vM) {
    if(b_==0) return 4.450000;
    if(b_==1) return 3.350000;
    if(b_==2) return 1.650000;
  };

  if(v_==vPt) {
    if(b_==0) return 1.850000;
    if(b_==1) return 4.450000;
  };


  fprintf(stderr,"ERROR: Binning::GetAziMax needs to be updated!\n");
  return 5;
};

Binning::~Binning() {};

