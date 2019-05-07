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
#include "TCanvas.h"
#include "TLorentzVector.h"
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
    void CalculateAsymmetries();
    void EvalAsymmetry(TGraphErrors * asymGr, TH1D * mdistL, TH1D * mdistR);

    void DrawBoundLines();
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
    static const Int_t nBinsMax = 7;
    Float_t minIV[nIV];
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
    TH1D * wDist1[nIV][nBinsMax]; // fill appropriate IV (finely binned)
    TH1D * mDist1[nSpin][nIV][nBinsMax]; // fill, e.g., Sin(phiR) for each spin
    TGraphErrors * asym1[nIV][nBinsMax]; // asymmetry
    TH1D * bDist1[nIV]; // full IV distribution, for showing bin bounds
    TCanvas * bDistCanv1[nIV]; // canvas for bDist


    // 2-dim binning:
    // distributions for each pair of bins in 2 IVs; the third
    // IV not considered is integrated over
    // [spin] [IV1] [IV2] [bin1] [bin2]
    // -- [IV1] is the first IV
    // -- [IV2] is the second IV
    // -- [bin1] is the first IV's bin
    // -- [bin2] is the second IV's bin
    TH2D * wDist2[nIV][nIV][nBinsMax][nBinsMax];
    TH1D * mDist2[nSpin][nIV][nIV][nBinsMax][nBinsMax];
    TGraphErrors * asym2[nIV][nIV][nBinsMax][nBinsMax];
    TH2D * bDist2[nIV][nIV];
    TCanvas * bDistCanv2[nIV][nIV];

    // 3-dim binning:
    // distributions for each triple of bins; since no IVs are integated
    // over, we only need one distribution; the three bins are in alphabetical
    // order
    // [vM bin] [vX bin] [vZ bin]
    TH3D * wDist3[nBinsMax][nBinsMax][nBinsMax];
    TH1D * mDist3[nSpin][nBinsMax][nBinsMax][nBinsMax];
    TGraphErrors * asym3[nBinsMax][nBinsMax][nBinsMax];


    // bin boundaries
    Int_t nBins[nIV]; // the number of bins
    Float_t bound[nIV][nBinsMax]; // [IV] [bin#]; this is the lower bound
                                  // for bin "bin#"; to get the upper bound
                                  // of the highest bin, use bin#+1
    Float_t nModMax;

    TLine * boundLine1[nIV][nBinsMax];
    TLine * vertLine[nIV][nIV][nBinsMax];
    TLine * horizLine[nIV][nIV][nBinsMax];

    TString IVname[nIV];
    TString IVtitle[nIV];
    TString ModulationTitle;
    
  private:
    TString SpinName[nSpin];
    TString SpinTitle[nSpin];
    TString plotTitle,plotName;
    TString canvName;


    Float_t angle;
    Float_t iv[nIV]; // independent variable
    Float_t modulation;
    Int_t binn[nIV];
    Int_t spinn;
    Int_t pointCnt;

    Double_t yL,yR;
    Double_t asymNumer,asymDenom;
    Double_t asymVal,modVal;
    Double_t asymErr,modErr;

    Double_t bMax;



  ClassDef(Asymmetry,1);
};

#endif
