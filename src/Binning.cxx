#include "Binning.h"

ClassImp(Binning)

using namespace std;


Binning::Binning(Int_t pairType_) {
  //printf("Instantiating Binning...\n");

  // get hadron indices (bin bounds depends on hadron type)
  DecodePairType(pairType_,whichHad[qA],whichHad[qB]);
  numKaons = 0;
  for(int h=0; h<2; h++) { if(whichHad[h]==kKp || whichHad[h]==kKm) numKaons++; };

  // set minimum and maximum IV values
  minIV[vM] = 0;   maxIV[vM] = 3;
  minIV[vX] = 0;   maxIV[vX] = 1;
  minIV[vZ] = 0;   maxIV[vZ] = 1;
  minIV[vPt] = 0;  maxIV[vPt] = 2;
  minIV[vPh] = 0;  maxIV[vPh] = 10;
  minIV[vQ] = 0;   maxIV[vQ] = 12;
  for(int v=0; v<nIV; v++) nBins[v]=-1;


  // set minimum bin boundaries
  for(int v=0; v<nIV; v++) AddBinBound(v,minIV[v]);


  // set main bin boundaries
  if(numKaons==0) {

    // -- M_h (dihadron invariant mass)
    AddBinBound(vM,0.53); // 5 quantiles (from GetQuantiles.C)
    AddBinBound(vM,0.71);
    AddBinBound(vM,0.84);
    AddBinBound(vM,1.03);

    // -- x (bjorken-x)
    AddBinBound(vX,0.22); // 5 quantiles (from GetQuantiles.C)
    AddBinBound(vX,0.26);
    AddBinBound(vX,0.32);
    AddBinBound(vX,0.39);

    // -- z (fragmentation fraction)
    AddBinBound(vZ,0.47); // 5 quantiles (from GetQuantiles.C)
    AddBinBound(vZ,0.55);
    AddBinBound(vZ,0.63);
    AddBinBound(vZ,0.73);

    // -- PhPerp (transverse momentum of dihadron)
    AddBinBound(vPt,0.29); // 5 quantiles (from GetQuantiles.C)
    AddBinBound(vPt,0.42);
    AddBinBound(vPt,0.55);
    AddBinBound(vPt,0.70);

    // -- Ph (magnitude of momentum sum of dihadron)
    AddBinBound(vPh,3.15); // 5 quantiles (from GetQuantiles.C)
    AddBinBound(vPh,3.60);
    AddBinBound(vPh,4.10);
    AddBinBound(vPh,4.80);

    // -- Q^2
    AddBinBound(vQ,2.84); // 5 quantiles (from GetQuantiles.C)
    AddBinBound(vQ,3.33);
    AddBinBound(vQ,3.92);
    AddBinBound(vQ,4.85);
    
  } else if(numKaons==1) {

    // -- mass
    AddBinBound(vM,0.85);
    AddBinBound(vM,1.1);

    // -- x
    AddBinBound(vX,0.2);
    AddBinBound(vX,0.3);

    // -- z
    AddBinBound(vZ,0.5);
    AddBinBound(vZ,0.7);

    // -- other IV binning schemes can be added later when kaons are included
  };

  // set maximum bin boundaries
  for(int v=0; v<nIV; v++) AddBinBound(v,maxIV[v]);


  // set IV names and titles
  IVname[vM] = "M";
  IVname[vX] = "X";
  IVname[vZ] = "Z";
  IVname[vPt] = "Pt";
  IVname[vPh] = "Ph";
  IVname[vQ] = "Q";

  IVtitle[vM] = "M_{h}";
  IVtitle[vX] = "x";
  IVtitle[vZ] = "z";
  IVtitle[vPt] = "P_{h}^{perp}";
  IVtitle[vPh] = "P_{h}";
  IVtitle[vQ] = "Q^{2}";


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


Binning::~Binning() {};
