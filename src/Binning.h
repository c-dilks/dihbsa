#ifndef Binning_
#define Binning_

#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <map>
#include <vector>

// ROOT
#include "TSystem.h"
#include "TObject.h"
#include "TTree.h"
#include "TChain.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TString.h"
#include "TMath.h"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TGraphErrors.h"
#include "TLine.h"

// dihbsa
#include "Constants.h"
#include "Tools.h"


class Binning : public TObject
{
  public:
    Binning() {}; // empty default constructor for streaming instances to ROOT files
    Binning(Int_t pairType_);
    ~Binning();
    void AddBinBound(Int_t ivIdx, Float_t newBound);
    void PrintBinBounds();
    Int_t GetBin(Int_t v_, Float_t iv_);
    TString GetBoundStr(Int_t v_, Int_t b_);

    Bool_t ValidIV(Int_t v_) { return v_>=0 && v_<nIV; };



    // enumerators 
    enum ivEnum { vX, vM, vZ, vPt, vPh, vQ, nIV }; // Independent Variables (IV)
    Float_t minIV[nIV];
    Float_t maxIV[nIV];


    // bin boundaries
    Int_t nBins[nIV]; // the number of bins
    std::vector<Float_t> bound[nIV]; // this is the lower bound
                                // for bin "bin#"; to get the upper bound
                                // of the highest bin, use bin#+1


    TString IVname[nIV];
    TString IVtitle[nIV];

    
  private:
    Int_t whichHad[2];
    Int_t numKaons;


  ClassDef(Binning,1);
};

#endif
