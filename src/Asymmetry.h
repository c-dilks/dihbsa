#ifndef Asymmetry_
#define Asymmetry_

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
#include "TLorentzVector.h"
#include "TString.h"
#include "TMath.h"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"

// dihbsa
#include "Constants.h"
#include "Tools.h"
#include "Trajectory.h"
#include "DIS.h"


class Asymmetry : public TObject
{
  public:
    Asymmetry(Int_t phiModulation, Bool_t singleBinMode);
    ~Asymmetry();
    void AddBinBound(Int_t ivIdx, Float_t newBound);
    void PrintBinBounds();
    Int_t GetBin(Int_t v_, Float_t iv_);

    void Write(TFile * f);

    void FillPlots();
    Float_t EvalModulation(Float_t PhiH_, Float_t PhiR_);
    Int_t SpinState(Int_t spin_);
    void ResetVars();

    Int_t nEvents;


    // modulations
    enum modEnum {
      modSinPhiR,
      modSinPhiHR,
      nMod
    };
    Int_t whichMod;

    // event-level variables -- these must be set for each event,
    // prior to calling FillHists
    Float_t Mh;
    Float_t x;
    Float_t z;
    Int_t eSpin;
    Int_t pSpin;
    Float_t PhiH;
    Float_t PhiR;
    Float_t PhPerp;


    // enumerators 
    enum ivEnum { vM, vX, vZ, nIV }; // Independent Variables (IV)
    enum spinEnum { sP, sM, nSpin }; // spin state
    static const Int_t nBinsMax = 20;
    Float_t maxIV[nIV];

    static const Int_t w1Bins = 100; // number of bins for wDist1 plots
    static const Int_t w2Bins = 50; // number of bins for wDist2 plots
    static const Int_t w3Bins = 30; // number of bins for wDist3 plots
    static const Int_t nModBins = 7; // number of bins in azimuthal modulation
    Float_t modMax;


    // 1-dim binning:
    // distributions for each bin for a specific IV; 
    // [spin] [IV] [bin]
    // -- [IV] is that specific IV; the other 2 IVs are integrated 
    // -- [bin] is the bin in that IV
    TH1D * mDist1[nSpin][nIV][nBinsMax]; // filled with, e.g., Sin(phiR) for each spin
    TH1D * wDist1[nIV][nBinsMax]; // filled with appropriate IV (finely binned)


    // 2-dim binning:
    // distributions for each pair of bins in 2 IVs; the third
    // IV not considered is integrated over
    // [spin] [IV1] [IV2] [bin1] [bin2]
    // -- [IV1] is the first IV
    // -- [IV2] is the second IV
    // -- [bin1] is the first IV's bin
    // -- [bin2] is the second IV's bin
    TH2D * mDist2[nSpin][nIV][nIV][nBinsMax][nBinsMax];
    TH1D * wDist2[nIV][nIV][nBinsMax][nBinsMax];

    // 3-dim binning:
    // distributions for each triple of bins; since no IVs are integated
    // over, we only need one distribution; the three bins are in alphabetical
    // order
    // [vM bin] [vX bin] [vZ bin]
    TH3D * mDist3[nSpin][nBinsMax][nBinsMax][nBinsMax];
    TH1D * wDist3[nBinsMax][nBinsMax][nBinsMax];


    // bin boundaries
    Int_t nBins[nIV]; // the number of bins
    Float_t bound[nIV][nBinsMax]; // [IV] [bin#]; this is the lower bound
                                  // for bin "bin#"; to get the upper bound
                                  // of the highest bin, use bin#+1
    Float_t nModMax;


    
  private:
    TString IVname[nIV];
    TString IVtitle[nIV];
    TString SpinName[nSpin];
    TString SpinTitle[nSpin];
    TString plotTitle,plotName;
    TString ModulationTitle;


    Float_t angle;
    Float_t iv[nIV]; // independent variable
    Float_t modulation;
    Int_t binn[nIV];
    Int_t spinn;



  ClassDef(Asymmetry,1);
};

#endif
