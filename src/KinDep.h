#ifndef KinDep_
#define KinDep_

#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <sstream>
#include <stdexcept>

// ROOT
#include "TSystem.h"
#include "TObject.h"
#include "TFile.h"
#include "TString.h"
#include "TMath.h"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TF1.h"
#include "TGraphErrors.h"
#include "TCanvas.h"
#include "TLine.h"

// dihbsa
#include "Constants.h"
#include "Tools.h"
#include "Asymmetry.h"


class KinDep : public TObject
{
  public:
    KinDep(Asymmetry * asym_);
    ~KinDep();
    void FillAsymGraphs();
    void SetAsymPoint(
      TGraphErrors * asymGr_, TGraphErrors * asym_, TH1 * wdist_,
      Int_t &cnt_, Int_t wdistAxis
    );
    void FillCanvases();
    void FormatAsymGr(TGraphErrors * g, Int_t ivNum);
    void WriteObjects(TFile * f);
    void PrintPNGs(Int_t whichPhiR_);

    Asymmetry * A;
    Int_t canvSize;
    TString fitName;

    
    // 1D canvases
    // divided into three panels, to show asym
    // plotted against each IV
    TCanvas * canv1;
    TGraphErrors * asymGr1[Asymmetry::nIV];

    // 2D canvases
    // [iv0] [iv1] [iv1bin]
    // line of panels, each showing asym vs. [iv0]
    // for a bin in [iv1]
    TCanvas * canv2[Asymmetry::nIV][Asymmetry::nIV];
    TGraphErrors * asymGr2[Asymmetry::nIV][Asymmetry::nIV][Asymmetry::nBinsMax];

    // 3D canvases
    // [iv0] [iv1] [iv2] [iv1bin] [iv2bin]
    // grid of panels, each showing asym vs. [iv0]
    // for a bin in ([iv1],[iv2])
    TCanvas * canv3[Asymmetry::nIV][Asymmetry::nIV][Asymmetry::nIV];
    TGraphErrors * asymGr3[Asymmetry::nIV][Asymmetry::nIV][Asymmetry::nIV][Asymmetry::nBinsMax][Asymmetry::nBinsMax];


    TLine * zeroLine[Asymmetry::nIV];
    Float_t ivMin[Asymmetry::nIV];
    Float_t ivMax[Asymmetry::nIV];


    
  private:
    Int_t N;
    Int_t NB[Asymmetry::nIV];
    TString canvName;
    TString grName,grTitle;

    TF1 * fitFunc;
    Float_t asymValue,asymError;
    Float_t kinValue,kinError;

    TGraphErrors * asymCurr;
    TH3D * wDistCurr;

    TString printName;
    TString printSuffix;

    Int_t cnt;

  ClassDef(KinDep,1);
};

#endif
