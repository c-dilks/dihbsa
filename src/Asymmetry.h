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
#include "TGraphErrors.h"
#include "TLine.h"

// dihbsa
#include "Constants.h"
#include "Tools.h"
#include "Trajectory.h"
#include "DIS.h"
#include "Binning.h"


class Asymmetry : public TObject
{
  public:
    Asymmetry(
      Binning * binScheme,
      Int_t phiModulation, Int_t dimension, 
      Int_t var0=0,  Int_t bin0=0,
      Int_t var1=-1, Int_t bin1=-1,
      Int_t var2=-1, Int_t bin2=-1
    );
    ~Asymmetry();

    void CalculateAsymmetries();
    void EvalAsymmetry(TGraphErrors * asymGr, TH1D * mdistL, TH1D * mdistR);


    void FillPlots();
    Float_t EvalModulation();
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
    Int_t whichDim;
    Binning * BS;


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
    Float_t pol;


    // number of bins
    static const Int_t w1Bins = 100; // number of bins for wDist1 plots
    static const Int_t w2Bins = 50; // number of bins for wDist2 plots
    static const Int_t w3Bins = 30; // number of bins for wDist3 plots
    static const Int_t nModBins = 7; // number of bins in azimuthal modulation
    Float_t modMax;


    // "iv dist": finely-binned IV distribution (for each whichDim)
    TH1D * ivDist1;
    TH2D * ivDist2;
    TH3D * ivDist3;
    TString ivName,ivTitle;

    // "azimuthal modulation dist" filled with, e.g., Sin(phiR) for each spin, binned
    // for the asymmetry plots
    TH1D * aziDist[nSpin];
    TString aziName[nSpin];
    TString aziTitle[nSpin];

    // "finely-binnd azimuthal modulation dist" filled for all spins, used for getting
    // mean value of azimuthal modulation 
    TH1D * modBinDist[nModBins]; // one for each modulation bin
    TH1D * modDist; // for all modulation bins
    TH2D * IVvsModDist;
    TString modName,modTitle;
    TString modBinName[nModBins];
    TString modBinTitle[nModBins];
    TString IVvsModName,IVvsModTitle;

    // asymmetry vs. azimuthal modulation bin
    TGraphErrors * asymGr;
    TString asymName,asymTitle;
    
    // bin lines
    TLine * boundLine;

    TString ModulationTitle;

    Bool_t debug;
    Bool_t success;
    Bool_t successIVmode;
    
  private:
    // variables for each dimension
    Int_t I[3]; // IV index number
    Int_t B[3]; // IV bin number
    Float_t iv[3]; // IV
    Float_t ivMin[3]; // IV minimum value (for histo ranges)
    Float_t ivMax[3]; // IV maximum value (for histo ranges)
    TString ivN[3]; // IV name
    TString ivT[3]; // IV title
    TString binT,binN; // bin title/name suffixes

    Float_t modulation;
    Int_t modbin;
    Int_t spinn;
    Int_t pointCnt;

    Double_t rellumNumer,rellumDenom,rellum;

    Double_t yL,yR;
    Double_t asymNumer,asymDenom;
    Double_t asymVal,modVal;
    Double_t asymErr,modErr;

    Double_t bMax;




  ClassDef(Asymmetry,1);
};

#endif
