// asymmetry fitter
// - fits spinroot cat file (from catSpinroot.cpp) for extracting asymmetry amplitudes
// - produces the file spinroot/asym.root
// - you can specify a specific `spinroot` directory, if you want

#include <cstdlib>
#include <iostream>
#include <vector>
#include <map>

// ROOT
#include "TFile.h"
#include "TString.h"
#include "TMath.h"
#include "TRegexp.h"
#include "TH1.h"
#include "TH2.h"
#include "TF1.h"
#include "TCanvas.h"
#include "TLine.h"
#include "TStyle.h"
#include "TMultiGraph.h"
#include "TSystem.h"

// DihBsa
#include "Constants.h"
#include "Binning.h"
#include "Asymmetry.h"


// global variables
Int_t N_AMP,N_D;



int main(int argc, char** argv) {

  // ARGUMENTS
  TString spinrootDir = "spinroot";
  if(argc>1) spinrootDir = TString(argv[1]);

  gStyle->SetOptFit(1);

  // open spinroot cat file and result file
  TFile * catFile = new TFile(TString(spinrootDir+"/cat.root"),"READ");
  TFile * asymFile = new TFile(TString(spinrootDir+"/asym.root"),"RECREATE");


  // instantiate Binning and Asymmetry, and extract them from catFile
  Binning * BS;
  Asymmetry * A;
  std::map<Int_t, Asymmetry*> asymMap;
  catFile->GetObject("BS",BS);
  for(Int_t bn : BS->binVec) {
    A = new Asymmetry(BS,bn);
    if(A->success) {
      A->AppendData(catFile);
      asymMap.insert(std::pair<Int_t,Asymmetry*>(bn,A));
    }
    else return 0;
  };
  TString dihTitle = PairTitle(BS->whichHad[qA],BS->whichHad[qB]);
  TString dihName = PairName(BS->whichHad[qA],BS->whichHad[qB]);


  // build maps from Binning::binVec number to plots etc.
  TGraphErrors * kindepGr;
  TGraphErrors * RFkindepGr[Asymmetry::nAmp];
  TGraphErrors * chindfGr;
  TGraphErrors * rellumGr;
  TMultiGraph * multiGr;

  std::map<Int_t, TGraphErrors*> kindepMap;
  std::map<Int_t, TGraphErrors*> RFkindepMap[Asymmetry::nAmp];
  std::map<Int_t, TMultiGraph*> multiMap;
  std::map<Int_t, TGraphErrors*> chindfMap;
  std::map<Int_t, TGraphErrors*> rellumMap;

  TString grTitle,grTitleSuffix,grName,grNameSuffix;
  Bool_t first = true;
  for(Int_t bn : BS->binVec) {
    A = asymMap.at(bn);

    // get number of RooFit params
    if(first) { 
      N_AMP = A->nAmpUsed;
      N_D = A->nDparamUsed;
      first = false;
    };

    // set graph title and name suffixes
    switch(BS->dimensions) {
      case 1:
        grTitleSuffix = BS->GetIVtitle(0) + ";" + BS->GetIVtitle(0);
        grNameSuffix = BS->GetIVname(0);
        break;
      case 2:
        grTitleSuffix = BS->GetIVtitle(0) + " :: " +
          BS->GetBoundStr(bn,1) + ";" +
          BS->GetIVtitle(0);
        grNameSuffix = Form("%s_bin_%s%d",
            (BS->GetIVname(0)).Data(),
            (BS->GetIVname(1)).Data(), BS->UnhashBinNum(bn,1)
            );
        break;
      case 3:
        grTitleSuffix = BS->GetIVtitle(0) + " :: " +
          BS->GetBoundStr(bn,1) + ", " +
          BS->GetBoundStr(bn,2) + ";" +
          BS->GetIVtitle(0);
        grNameSuffix = Form("%s_bin_%s%d_%s%d",
            (BS->GetIVname(0)).Data(),
            (BS->GetIVname(1)).Data(), BS->UnhashBinNum(bn,1),
            (BS->GetIVname(2)).Data(), BS->UnhashBinNum(bn,2)
            );
        break;
    };


    // instantiate graphs; only needs to be done for each IV1 and IV2 bin, since the
    // horizontal axis of these graphs are all IV0 (hence the if statement here)
    if(BS->UnhashBinNum(bn,0)==0) {

      // instantiate kindep graph, for linear fit ("l.f.") result
      grTitle = dihTitle + " A_{LU}[" + A->ModulationTitle + "]_{l.f.} " + 
        " vs. " + grTitleSuffix;
      grName = "kindep_" + grNameSuffix;
      kindepGr = new TGraphErrors();
      kindepGr->SetName(grName);
      kindepGr->SetTitle(grTitle);

      // instantiate multiGraph, for plotting kindep graphs together
      grTitle = dihTitle + " A_{LU}[" + A->ModulationTitle + "] " + 
        " vs. " + grTitleSuffix;
      grName = "multiGr_" + grNameSuffix;
      multiGr = new TMultiGraph();
      multiGr->SetName(grName);
      multiGr->SetTitle(grTitle);

      // instantiate kindep graphs for maximum likelihood method (m.l.m.) result
      // for each fit parameter
      for(int aa=0; aa<N_AMP; aa++) {
        grTitle = dihTitle + " " + TString(A->rfA[aa]->GetTitle()) + "_{m.l.m.} " +
          " vs. " + grTitleSuffix;
        grName = "RF_A" + TString::Itoa(aa,10) + "_kindep_" + grNameSuffix;
        RFkindepGr[aa] = new TGraphErrors();
        RFkindepGr[aa]->SetName(grName);
        RFkindepGr[aa]->SetTitle(grTitle);
      };

      // instantiate chi2/ndf graphs
      grTitle = "#chi^{2}/NDF of " +
        dihTitle + " A_{LU}[" + A->ModulationTitle + "]_{l.f.} " + 
        " vs. " + grTitleSuffix;
      grName = "chindf_" + grNameSuffix;
      chindfGr = new TGraphErrors();
      chindfGr->SetName(grName);
      chindfGr->SetTitle(grTitle);

      // instantiate relative luminosity graphs
      grTitle = "relative luminosity vs. " + grTitleSuffix;
      grName = "rellum_" + grNameSuffix;
      rellumGr = new TGraphErrors();
      rellumGr->SetName(grName);
      rellumGr->SetTitle(grTitle);

    }


    // insert objects into maps (note: these are many-to-one maps, i.e., several
    // bin numbers will map to the same pointer)
    kindepMap.insert(std::pair<Int_t,TGraphErrors*>(bn,kindepGr));
    multiMap.insert(std::pair<Int_t,TMultiGraph*>(bn,multiGr));
    chindfMap.insert(std::pair<Int_t,TGraphErrors*>(bn,chindfGr));
    rellumMap.insert(std::pair<Int_t,TGraphErrors*>(bn,rellumGr));
    for(int aa=0; aa<N_AMP; aa++) {
      RFkindepMap[aa].insert(std::pair<Int_t,TGraphErrors*>(bn,RFkindepGr[aa]) );
    };
  }; // eo binVec loop
};

