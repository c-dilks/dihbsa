#ifndef Asymmetry_
#define Asymmetry_

#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <math.h>
#include <map>
#include <vector>
#include <vector>

// ROOT
#include "TSystem.h"
#include "TStyle.h"
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
#include "TF1.h"
#include "TF2.h"
#include "TGraphErrors.h"
#include "TGraph2DErrors.h"
#include "TLine.h"

// RooFit
#include <RooGlobalFunc.h>
#include <RooGenericPdf.h>
#include <RooFitResult.h>
#include <RooExtendPdf.h>
#include <RooAbsReal.h>
#include <RooArgSet.h>
#include <RooDataSet.h>
#include <RooRealVar.h>
#include <RooCategory.h>
#include <RooDataHist.h>
#include <RooSimultaneous.h>
#include <RooAddPdf.h>
#include <RooPlot.h>
#include <RooNLLVar.h>


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
      Int_t phiModulation, Int_t dimension=1, 
      Int_t var0=0,  Int_t bin0=0,
      Int_t var1=-1, Int_t bin1=-1,
      Int_t var2=-1, Int_t bin2=-1
    );
    ~Asymmetry();

    Bool_t InitRooFit();
    void CalculateAsymmetries();
    void CalculateRooAsymmetries();
    void SetAsymGrPoint(Int_t modBin_, Int_t modBin2_=-1);


    Bool_t FillPlots();
    Float_t EvalModulation();
    Float_t EvalWeight();
    Int_t SpinState(Int_t spin_);

    void ResetVars();
    void PrintSettings();

    Double_t nEvents;
    Double_t yield[nSpin];

    Bool_t success;
    Bool_t successIVmode;
    Bool_t debug;


    // modulations
    enum modEnum {
      modSinPhiR,
      modSinPhiHR,
      weightSinPhiHR,
      modSinPhiH,
      mod2d,
      nMod
    };
    Int_t whichMod;
    Int_t whichDim;
    Bool_t asym2d;
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
    Float_t theta;
    Float_t pol;


    // number of bins
    static const Int_t iv1Bins = 100; // number of bins for ivDist1 plots
    static const Int_t iv2Bins = 50; // number of bins for ivDist2 plots
    static const Int_t iv3Bins = 30; // number of bins for ivDist3 plots
    static const Int_t nModBins = 7; // number of bins in azimuthal modulation
    static const Int_t nModBins2 = 5; // number of bins in 2d azimuthal modulation

    Float_t modMaxDefault;
    Float_t modMax,aziMax;
    Float_t weight;


    // "iv dist": finely-binned IV distribution (for each whichDim)
    TH1D * ivDist1;
    TH2D * ivDist2;
    TH3D * ivDist3;
    TString ivName,ivTitle;

    // "azimuthal modulation dist" filled with, e.g., Sin(phiR) for each spin, binned
    // for the asymmetry plots
    TH1D * aziDist[nSpin];
    TH2D * aziDist2[nSpin]; // for 2d modulations
    TString aziName[nSpin];
    TString aziTitle[nSpin];

    // "finely-binnd azimuthal modulation dist" filled for all spins, used for getting
    // mean value of azimuthal modulation 
    TH1D * modBinDist[nModBins]; // one for each modulation bin
    TH1D * modDist; // for all modulation bins
    TH2D * IVvsModDist;
    TString modName,modTitle;
    TString modBinName;
    TString modBinTitle;
    TString IVvsModName,IVvsModTitle;
    // for 2d modulations:
    TH2D * modBinDist2[nModBins2][nModBins2];
    TH2D * modDist2;

    // asymmetry vs. azimuthal modulation bin
    TGraphErrors * asymGr;
    TString asymName,asymTitle;
    TF1 * fitFunc;
    TString fitFuncName;
    // for 2d modulations:
    TGraph2DErrors * asymGr2;
    TH2D * asymGr2hist;
    TF2 * fitFunc2;

    // bin lines
    TLine * boundLine;

    TString ModulationTitle,ModulationName;



    
    // variables for each dimension
    Int_t I[3]; // IV index number
    Int_t B[3]; // IV bin number
    Float_t iv[3]; // IV
    Float_t ivMin[3]; // IV minimum value (for histo ranges)
    Float_t ivMax[3]; // IV maximum value (for histo ranges)
    TString ivN[3]; // IV name
    TString ivT[3]; // IV title
    TString binT,binN; // bin title/name suffixes

    Double_t rNumer,rDenom,rellum,rellumErr;



    // RooFit variables
    Bool_t roofitter;
    RooDataSet * rfData;
    RooGenericPdf * rfPdf[2];
    //RooExtendPdf * rfModelExt[2];
    //RooAddPdf * rfPdf[2];
    RooSimultaneous * rfSimPdf;
    RooCategory * rfSpinCateg;
    RooFitResult * rfResult;
    RooArgSet * rfVars;
    RooArgSet * rfParams[2];
    TString rfPdfFormu[2];
    TString rfSpinName[2];

    TString rfModulation[nMod];
    TString pwFactorSP,pwFactorPP;
    TString asymExpansion;
    TString preFactor[2];
    Float_t rfParamRange;


    RooRealVar *rfPhiH, *rfPhiR, *rfTheta;
    RooRealVar *rfWeight;
    static const Int_t nAmp = 4;
    Int_t nAmpUsed;
    TString rfAname[nAmp];
    RooRealVar *rfA[nAmp];
    RooRealVar *rfYieldBoth;
    RooRealVar *rfYield[2];

    RooNLLVar * rfNLL;
    RooPlot * rfNLLplot[nAmp];

  private:

    Float_t modulation;
    Int_t modbin,modbinH,modbinR;
    Int_t spinn;
    Int_t pointCnt;


    Double_t yL,yR;
    Double_t asymNumer,asymDenom;
    Double_t asymVal,modVal,modValH,modValR;
    Double_t asymErr,modErr,modErrH,modErrR;

    Double_t bMax;






  ClassDef(Asymmetry,1);
};

#endif
