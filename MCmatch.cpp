#include <cstdlib>
#include <iostream>

// ROOT
#include "TFile.h"
#include "TString.h"
#include "TMath.h"
#include "TRegexp.h"
#include "TH1.h"
#include "TH2.h"
#include "TProfile.h"
#include "TCanvas.h"

// DihBsa
#include "Constants.h"
#include "Tools.h"
#include "DIS.h"
#include "Dihadron.h"
#include "Trajectory.h"
#include "EventTree.h"

Float_t Delta(Float_t vGen, Float_t vRec, Bool_t adjAngle=false);

int main(int argc, char** argv) {

  // ARGUMENTS
  TString inDir = "outroot.MC.rec";
  Int_t whichPair = EncodePairType(kPip,kPim);
  if(argc>1) inDir = TString(argv[1]);
  if(argc>2) whichPair = (Int_t)strtof(argv[2],NULL);


  // open input files and define output file
  EventTree * ev = new EventTree(TString(inDir+"/*.root"),whichPair);
  TFile * outfile = new TFile("match.root","RECREATE");

  Int_t whichHad[2];
  TString hadName[2];
  TString hadTitle[2];
  int h;
  DecodePairType(whichPair,whichHad[qA],whichHad[qB]);
  for(int h=0; h<2; h++) {
    hadName[h] = PairHadName(whichHad[qA],whichHad[qB],h);
    hadTitle[h] = PairHadTitle(whichHad[qA],whichHad[qB],h);
  };


  // define histograms
  const Int_t NBINS = 100;
  const Float_t Dlim = 0.3;
  TH1F * Ddist = new TH1F("Ddist","D distribution",NBINS,0,10);
  TH1F * DdistZoom = new TH1F("DdistZoom","D distribution (zoom)",NBINS,0,Dlim);

  TH2F * hadEDiffVsD[2];
  TH2F * hadPtDiffVsD[2];
  TH2F * hadPhiDiffVsD[2];
  TH2F * hadThetaDiffVsD[2];
  for(h=0; h<2; h++) {
    hadEDiffVsD[h] = new TH2F(
      TString(hadName[h]+"EDiffVsD"),
      TString(hadTitle[h]+" #Delta E vs. D"),
      NBINS,0,Dlim,NBINS,-1,1);
    hadPtDiffVsD[h] = new TH2F(
      TString(hadName[h]+"PtDiffVsD"),
      TString(hadTitle[h]+" #Delta p_{T} vs. D"),
      NBINS,0,Dlim,NBINS,-1,1);
    hadPhiDiffVsD[h] = new TH2F(
      TString(hadName[h]+"PhiDiffVsD"),
      TString(hadTitle[h]+" #Delta #phi vs. D"),
      NBINS,0,Dlim,NBINS,-1,1);
    hadThetaDiffVsD[h] = new TH2F(
      TString(hadName[h]+"ThetaDiffVsD"),
      TString(hadTitle[h]+" #Delta #theta vs. D"),
      NBINS,0,Dlim,NBINS,-1,1);
  };


  // event loop
  printf("begin loop through %lld events...\n",ev->ENT);
  for(int i=0; i<ev->ENT; i++) {

    ev->GetEvent(i);

    if(ev->Valid()) {

      Ddist->Fill(ev->matchDiff);
      DdistZoom->Fill(ev->matchDiff);
      for(h=0; h<2; h++) {
        hadEDiffVsD[h]->Fill(ev->matchDiff,Delta(ev->gen_hadE[h],ev->hadE[h]));
        hadPtDiffVsD[h]->Fill(ev->matchDiff,Delta(ev->gen_hadPt[h],ev->hadPt[h]));
        hadPhiDiffVsD[h]->Fill(ev->matchDiff,Delta(ev->gen_hadPhi[h],ev->hadPhi[h]),1);
        hadThetaDiffVsD[h]->Fill(ev->matchDiff,Delta(ev->gen_hadTheta[h],ev->hadTheta[h]));
      };
    };
  };

  //write output
  Ddist->Write();
  DdistZoom->Write();
  for(h=0; h<2; h++) {
    hadEDiffVsD[h]->Write();
    hadPtDiffVsD[h]->Write();
    hadPhiDiffVsD[h]->Write();
    hadThetaDiffVsD[h]->Write();
  };

  outfile->Close();
};

Float_t Delta(Float_t vGen, Float_t vRec, Bool_t adjAngle) {
  if(vRec==0) {
    fprintf(stderr,"ERROR: vRec==0\n");
    return UNDEF;
  }
  if(adjAngle) return Tools::AdjAngle(vGen-vRec)/vRec;
  else return (vGen-vRec)/vRec;
};
