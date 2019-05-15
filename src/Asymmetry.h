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
#include "Binning.h"


class Asymmetry : public TObject
{
  public:
    Asymmetry(
      Binning * binScheme,
      Int_t phiModulation, Int_t dimension, 
      Int_t var0=0, Int_t var1=-1, Int_t var2=-1
    );
    ~Asymmetry();

    void CalculateAsymmetries();
    void EvalAsymmetry(TGraphErrors * asymGr, TH1D * mdistL, TH1D * mdistR);

    void DrawBoundLines();
    void WriteObjects(TFile * f);

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


    // enumerators 
    enum spinEnum { sP, sM, nSpin }; // spin state

    static const Int_t w1Bins = 100; // number of bins for wDist1 plots
    static const Int_t w2Bins = 50; // number of bins for wDist2 plots
    static const Int_t w3Bins = 30; // number of bins for wDist3 plots
    static const Int_t nModBins = 7; // number of bins in azimuthal modulation
    Float_t modMax;


    // "iv dist": finely-binned IV distribution (for each whichDim)
    vector<TH1D*> ivVec1;
    vector<TH2D*> ivVec2;
    vector<TH3D*> ivVec3;
    TH1D * ivDist1;
    TH2D * ivDist2;
    TH3D * ivDist3;
    TH1D * ivFullDist1;
    TH2D * ivFullDist2;
    TH3D * ivFullDist3;
    TCanvas * ivFullCanv;

    // "azimuthal modulation dist" filled with, e.g., Sin(phiR) for each spin
    vector<TH1D*> aziVec[nSpin]; 
    TH1D * aziDist[nSpin];

    // "finely-binnd azimuthal modulation dist" filled for all spins
    vector<TH1D*> modVec;
    TH1D * modDist;
    TH1D * modFullDist;
    TH2D * IVvsModDist;

    // asymmetry vs. azimuthal modulation bin
    vector<TGraphErrors*> asymVec;
    TGraphErrors * asymDist;
    
    // bin lines
    vector<TLine*> boundLineVec;
    TLine * boundLine;

    TString ModulationTitle;

    Bool_t success;
    Bool_t successIVmode;
    
  private:
    Int_t I[3];
    Float_t iv[3]; // independent variable
    Float_t ivMin[3];
    Float_t ivMax[3];
    Int_t nBin[3];
    TString ivN[3];
    TString ivT[3];

    TString SpinName[nSpin];
    TString SpinTitle[nSpin];
    TString plotTitle,plotName;
    TString canvName;


    Float_t angle;
    Float_t modulation;
    Int_t binn[nIV];
    Int_t spinn;
    Int_t pointCnt;

    Double_t rellumNumer,rellumDenom,rellum;

    Double_t yL,yR;
    Double_t asymNumer,asymDenom;
    Double_t asymVal,modVal;
    Double_t asymErr,modErr;

    Double_t bMax;

    Int_t fbin;



  ClassDef(Asymmetry,1);
};

#endif
