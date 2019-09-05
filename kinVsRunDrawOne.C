R__LOAD_LIBRARY(DihBsa)

#include "Constants.h"
#include "EventTree.h"
#include "Tools.h"
#include "TString.h"
#include "TMath.h"
#include "TFile.h"
#include "TCanvas.h"

void kinVsRunDrawOne(Int_t runnum=5032, TString infileN="vsRun.main.root") {

  //gStyle->SetOptStat(1101); // do not show num. entries
  gStyle->SetOptStat(0); // do not show num. entries

  TFile * infile = new TFile(infileN,"READ");
  TCanvas * canv1 = new TCanvas("canv1","canv1",1000,1000); canv1->Divide(4,2);
  TCanvas * canv2 = new TCanvas("canv2","canv2",1000,1000); canv2->Divide(3,1);
  
  TString runstr = Form("_run%d",runnum);

  TH1D * elePhiDist = (TH1D*) infile->Get(TString("RunByRun/elePhi/elePhi"+runstr));
  TH1D * xDist = (TH1D*) infile->Get(TString("RunByRun/X/X"+runstr));
  TH1D * MhDist = (TH1D*) infile->Get(TString("RunByRun/Mh/Mh"+runstr));
  TH1D * thetaDist = (TH1D*) infile->Get(TString("RunByRun/theta/theta"+runstr));
  TH1D * ZpairDist = (TH1D*) infile->Get(TString("RunByRun/Zpair/Zpair"+runstr));

  TH1D * hadEDist[2];
  TH1D * hadPtDist[2];
  TH1D * hadPhiDist[2];
  TString hadName[2];

  for(int h=0; h<2; h++) {
    hadName[h] = PairHadName(kPip,kPim,h);
    hadEDist[h] = (TH1D*) infile->Get(TString("RunByRun/"+hadName[h]+"hadE/"+hadName[h]+"hadE"+runstr));
    hadPtDist[h] = (TH1D*) infile->Get(TString("RunByRun/"+hadName[h]+"hadPt/"+hadName[h]+"hadPt"+runstr));
    hadPhiDist[h] = (TH1D*) infile->Get(TString("RunByRun/"+hadName[h]+"hadPhi/"+hadName[h]+"hadPhi"+runstr));
  };

  canv1->cd(1); xDist->Draw();
  canv1->cd(2); MhDist->Draw();
  canv1->cd(5); ZpairDist->Draw();
  canv1->cd(6); thetaDist->Draw();

  canv1->cd(3); hadEDist[qA]->Draw();
  canv1->cd(4); hadPtDist[qA]->Draw();
  canv1->cd(7); hadEDist[qB]->Draw();
  canv1->cd(8); hadPtDist[qB]->Draw();

  canv2->cd(1); elePhiDist->Draw();
  canv2->cd(2); hadPhiDist[qA]->Draw();
  canv2->cd(3); hadPhiDist[qB]->Draw();

};
