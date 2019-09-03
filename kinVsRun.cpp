#include <cstdlib>
#include <iostream>

// ROOT
#include "TFile.h"
#include "TTree.h"
#include "TChain.h"
#include "TString.h"
#include "TMath.h"
#include "TSystem.h"
#include "TRegexp.h"
#include "TH1.h"
#include "TH2.h"
#include "TCanvas.h"
#include "TGraph.h"

// DihBsa
#include "Constants.h"
#include "Tools.h"
#include "DIS.h"
#include "Trajectory.h"
#include "Dihadron.h"
#include "EventTree.h"

TString inDir;
Int_t whichPair;
Int_t whichHad[2];
TString hadName[2];
TString hadTitle[2];

void normalize(TH2D * d);

int main(int argc, char** argv) {

  // RUN NUMBER RANGE
  const Int_t runnumMin = 5000;
  const Int_t runnumMax = 5300;
  Int_t runnumBins = runnumMax - runnumMin;


  // ARGUMENTS
  inDir = "outroot.fall18";
  whichPair = EncodePairType(kPip,kPim);
  if(argc>1) inDir = TString(argv[1]);
  if(argc>2) whichPair = (Int_t)strtof(argv[2],NULL);


  // get hadron pair from whichPair; note that in the print out, the 
  // order of hadron 0 and 1 is set by Constants::dihHadIdx
  printf("whichPair = 0x%x\n",whichPair);
  DecodePairType(whichPair,whichHad[qA],whichHad[qB]);
  for(int h=0; h<2; h++) {
    hadName[h] = PairHadName(whichHad[qA],whichHad[qB],h);
    hadTitle[h] = PairHadTitle(whichHad[qA],whichHad[qB],h);
    printf("hadron %d:  idx=%d  name=%s  title=%s\n",
        h,dihHadIdx(whichHad[qA],whichHad[qB],h),hadName[h].Data(),hadTitle[h].Data());
  };

  EventTree * ev = new EventTree(TString(inDir+"/*.root"),whichPair);


  TFile * outfile = new TFile("vsRun.root","RECREATE");

  const Int_t NBINS = 100; // number of bins for kinematic variable
  Float_t deltaPhi;


  // DIS kinematics
  TH2D * WDist = new TH2D("WDist",
    "W distribution (w/o W cut) vs. runnum;runnum;W",
    runnumBins,runnumMin,runnumMax,NBINS,0,6);
  TH2D * XDist = new TH2D("XDist",
    "x distribution vs. runnum;runnum;x",
    runnumBins,runnumMin,runnumMax,NBINS,0,1);
  TH2D * YDist = new TH2D("YDist",
    "y distribution (w/o y cut) vs. runnum;runnum;y",
    runnumBins,runnumMin,runnumMax,NBINS,0,1);

  // electron kinematics
  TH2D * eleEDist = new TH2D("eleEDist",
    "e^{-} E distribution vs. runnum;runnum;E",
    runnumBins,runnumMin,runnumMax,NBINS,0,12);
  TH2D * elePtDist = new TH2D("elePtDist",
    "e^{-} p_{T} distribution vs. runnum;runnum;p_{T}",
    runnumBins,runnumMin,runnumMax,NBINS,0,4);
  TH2D * eleEtaDist = new TH2D("eleEtaDist",
    "e^{-} #eta distribution vs. runnum;runnum;#eta",
    runnumBins,runnumMin,runnumMax,NBINS,-3,6);
  TH2D * elePhiDist = new TH2D("elePhiDist",
    "e^{-} #phi distribution vs. runnum;runnum;#phi",
    runnumBins,runnumMin,runnumMax,NBINS,-PIe,PIe);


  // dihadron's hadron kinematics
  TH2D * hadEDist[2];
  TH2D * hadPDist[2];
  TH2D * hadPtDist[2];
  TH2D * hadEtaDist[2];
  TH2D * hadPhiDist[2];
  TH2D * hadZDist[2];
  for(int h=0; h<2; h++) {
    hadEDist[h] = new TH2D(TString(hadName[h]+"hadEDist"),
        TString(hadTitle[h]+" E vs. runnum;runnum;E"),
        runnumBins,runnumMin,runnumMax,NBINS,0,10);
    hadPDist[h] = new TH2D(TString(hadName[h]+"hadPDist"),
        TString(hadTitle[h]+" p vs. runnum;runnum;p"),
        runnumBins,runnumMin,runnumMax,NBINS,0,10);
    hadPtDist[h] = new TH2D(TString(hadName[h]+"hadPtDist"),
        TString(hadTitle[h]+" p_{T} vs. runnum;runnum;p_{T}"),
        runnumBins,runnumMin,runnumMax,NBINS,0,4);
    hadEtaDist[h] = new TH2D(TString(hadName[h]+"hadEtaDist"),
        TString(hadTitle[h]+" #eta vs. runnum;runnum;#eta"),
        runnumBins,runnumMin,runnumMax,NBINS,0,5);
    hadPhiDist[h] = new TH2D(TString(hadName[h]+"hadPhiDist"),
        TString(hadTitle[h]+" #phi vs. runnum;runnum;#phi"),
        runnumBins,runnumMin,runnumMax,NBINS,-PIe,PIe);
    hadZDist[h] = new TH2D(TString(hadName[h]+"hadZDist"),
        TString(hadTitle[h]+" z vs. runnum;runnum;z"),
        runnumBins,runnumMin,runnumMax,NBINS,0,1);
  };


  // dihadron kinematics
  TString plotTitle = "#Delta#phi = #phi(" + hadTitle[qA] + ")" +
    " - #phi(" + hadTitle[qB] + 
    ") distribution vs runnum;runnum;#Delta#phi";
  TH2D * deltaPhiDist = new TH2D("deltaPhiDist",
    plotTitle,
    runnumBins,runnumMin,runnumMax,NBINS,-PIe,PIe);

  TH2D * MhDist = new TH2D("MhDist",
    "M_{h} distribution vs. runnum;runnum;M_{h}",
    runnumBins,runnumMin,runnumMax,2*NBINS,0,3);
  TH2D * ZpairDist = new TH2D("ZpairDist",
    "z_{pair} distribution vs. runnum;runnum;z_{pair}",
    runnumBins,runnumMin,runnumMax,NBINS,0,1);
  TH2D * zetaDist = new TH2D("zetaDist",
    "#zeta distribution vs. runnum;runnum;#zeta",
    runnumBins,runnumMin,runnumMax,NBINS,-1,1);
  TH2D * xFDist = new TH2D("xFDist",
    "x_{F} distribution vs. runnum;runnum;x_{F}",
    runnumBins,runnumMin,runnumMax,NBINS,-2,2);
  TH2D * MmissDist = new TH2D("MmissDist",
    "M_{X} distribution vs. runnum;runnum;M_{X}",
    runnumBins,runnumMin,runnumMax,NBINS,-2,6);

  TH2D * PhiHDist = new TH2D("PhiHDist",
    "#phi_{h} distribution vs. runnum;runnum;#phi_{h}",
    runnumBins,runnumMin,runnumMax,NBINS,-PIe,PIe);
  TH2D * PhiRDist = new TH2D("PhiRDist",
    "#phi_{R} distribution vs. runnum;runnum;#phi_{R}",
    runnumBins,runnumMin,runnumMax,NBINS,-PIe,PIe);
  TH2D * PhiHRDist = new TH2D("PhiHRDist",
    "#phi_{h}-#phi_{R} distribution vs. runnum;runnum;#phi_{h}-#phi_{R}",
    runnumBins,runnumMin,runnumMax,NBINS,-PIe,PIe);


  // distributions for partial wave analysis
  TH2D * thetaDist = new TH2D("thetaDist",
    "#theta distribution vs. runnum;runnum;#theta",
    runnumBins,runnumMin,runnumMax,NBINS,0,PI);

  // multiplicities
  TH2D * partMultiplicity = new TH2D("partMultiplicity",
    "overall particle multiplicities (DIS cuts only) vs. runnum;runnum;particle",
    runnumBins,runnumMin,runnumMax,nParticles,0,nParticles);
  TH2D * obsMultiplicity = new TH2D("obsMultiplicity",
    "dihadrons' particle multiplicities vs. runnum;runnum;particle",
    runnumBins,runnumMin,runnumMax,nObservables,0,nObservables);
  for(int p=0; p<nParticles; p++) 
    partMultiplicity->GetYaxis()->SetBinLabel(p+1,PartTitle(p));
  for(int p=0; p<nObservables; p++) 
    obsMultiplicity->GetYaxis()->SetBinLabel(p+1,ObsTitle(p));


  // event-level distributions
  TH2D * helicityDist = new TH2D("helicityDist",
    "helicity vs. runnum;runnum;h",
    runnumBins,runnumMin,runnumMax,5,-2,3);




  // EVENT LOOP
  printf("begin loop through %lld events...\n",ev->ENT);
  Int_t hadI[2];
  for(int i=0; i<ev->ENT; i++) {

    ev->GetEvent(i);


    // fill multiplicity plots
    //------------------------
    // fill overall particle multiplicity
    if(ev->cutDIS) {
      for(int p=0; p<nParticles; p++) {
        if(ev->particleCnt[p]>0) partMultiplicity->Fill(ev->runnum,p,ev->particleCnt[p]);
      };
    };
    if(ev->cutDIS && ev->cutDihadronKinematics && ev->cutDiph[qA] && ev->cutDiph[qB]) {

      // fill observable multiplicity
      for(int h=0; h<2; h++) {
        hadI[h] = IO(ev->hadIdx[h]); // observable indices
        obsMultiplicity->Fill(ev->runnum,hadI[h]);
      };
    };


    // fill DIS kinematic plots
    // ------------------------
    if(ev->cutDihadron) {

      //ev->PrintEvent();

      if(ev->cutQ2 && ev->cutY) {
        WDist->Fill(ev->runnum,ev->W);
      };

      if(ev->cutQ2 && ev->cutW) {
        YDist->Fill(ev->runnum,ev->y);
      };
    };


    // fill dihadron kinematics plots
    // ------------------------------
    if(ev->Valid()) {

      eleEDist->Fill(ev->runnum,ev->eleE);
      elePtDist->Fill(ev->runnum,ev->elePt);
      eleEtaDist->Fill(ev->runnum,ev->eleEta);
      elePhiDist->Fill(ev->runnum,ev->elePhi);

      XDist->Fill(ev->runnum,ev->x);

      for(int h=0; h<2; h++) {
        hadEDist[h]->Fill(ev->runnum,ev->hadE[h]);
        hadPDist[h]->Fill(ev->runnum,ev->hadP[h]);
        hadPtDist[h]->Fill(ev->runnum,ev->hadPt[h]);
        hadEtaDist[h]->Fill(ev->runnum,ev->hadEta[h]);
        hadPhiDist[h]->Fill(ev->runnum,ev->hadPhi[h]);
        hadZDist[h]->Fill(ev->runnum,ev->Z[h]);
      };

      deltaPhi = Tools::AdjAngle(ev->hadPhi[qA] - ev->hadPhi[qB]);
      deltaPhiDist->Fill(ev->runnum,deltaPhi);

      MhDist->Fill(ev->runnum,ev->Mh);
      ZpairDist->Fill(ev->runnum,ev->Zpair);
      zetaDist->Fill(ev->runnum,ev->zeta);
      xFDist->Fill(ev->runnum,ev->xF);
      MmissDist->Fill(ev->runnum,ev->Mmiss);

      PhiHDist->Fill(ev->runnum,ev->PhiH);
      PhiRDist->Fill(ev->runnum,ev->PhiR);
      PhiHRDist->Fill(ev->runnum,ev->PhiHR);

      thetaDist->Fill(ev->runnum,ev->theta);

      helicityDist->Fill(ev->runnum,ev->helicity);
    };

  }; // eo event loop


  // normalize distributions
  normalize(WDist);
  normalize(XDist);
  normalize(YDist);

  normalize(eleEDist);
  normalize(elePtDist);
  normalize(eleEtaDist);
  normalize(elePhiDist);

  normalize(partMultiplicity);
  normalize(obsMultiplicity);

  for(int h=0; h<2; h++) {
    normalize(hadEDist[h]);
    normalize(hadPDist[h]);
    normalize(hadPtDist[h]);
    normalize(hadEtaDist[h]);
    normalize(hadPhiDist[h]);
    normalize(hadZDist[h]);
  };

  normalize(deltaPhiDist);
  normalize(MhDist);
  normalize(ZpairDist);
  normalize(zetaDist);
  normalize(xFDist);
  normalize(MmissDist);
  normalize(thetaDist);

  normalize(PhiHDist);
  normalize(PhiRDist);
  normalize(PhiHRDist);

  normalize(helicityDist);

  // write output
  WDist->Write();
  XDist->Write();
  YDist->Write();

  eleEDist->Write();
  elePtDist->Write();
  eleEtaDist->Write();
  elePhiDist->Write();

  partMultiplicity->Write();
  obsMultiplicity->Write();

  for(int h=0; h<2; h++) hadEDist[h]->Write();
  for(int h=0; h<2; h++) hadPDist[h]->Write();
  for(int h=0; h<2; h++) hadPtDist[h]->Write();
  for(int h=0; h<2; h++) hadEtaDist[h]->Write();
  for(int h=0; h<2; h++) hadPhiDist[h]->Write();
  for(int h=0; h<2; h++) hadZDist[h]->Write();

  deltaPhiDist->Write();
  MhDist->Write();
  ZpairDist->Write();
  zetaDist->Write();
  xFDist->Write();
  MmissDist->Write();
  thetaDist->Write();

  PhiHDist->Write();
  PhiRDist->Write();
  PhiHRDist->Write();


  helicityDist->Write();


  outfile->Close();
};


void normalize(TH2D * d) {
  printf("normalize %s\n",d->GetName());
  TH2D * c = (TH2D*)d->Clone();
  d->Reset();

  Int_t nx = c->GetNbinsX();
  Int_t ny = c->GetNbinsY();

  Double_t yl = c->GetYaxis()->GetXmin();
  Double_t yh = c->GetYaxis()->GetXmax();
  TH1D * pr = new TH1D("pr","pr",ny,yl,yh);

  Double_t bc,norm;
  for(int bx=1; bx<=nx; bx++) {
    
    pr->Reset();

    for(int by=1; by<=ny; by++) {
      //bn = c->GetBin(bx,by);
      //bc = c->GetBinContent(bn);
      bc = c->GetBinContent(bx,by);
      pr->SetBinContent(by,bc);
    };

    norm = pr->Integral();
    if(norm>0) pr->Scale(1/norm);

    for(int by=1; by<=ny; by++) {
      bc = pr->GetBinContent(by);
      d->SetBinContent(bx,by,bc);
    };
  };

  delete pr;
   
};
      

  
