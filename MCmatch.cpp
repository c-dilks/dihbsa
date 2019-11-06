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
#include "TProfile.h"
#include "TCanvas.h"

// DihBsa
#include "Constants.h"
#include "Tools.h"
#include "DIS.h"
#include "Dihadron.h"
#include "Trajectory.h"
#include "EventTree.h"

int main(int argc, char** argv) {

   // ARGUMENTS
   TString inDir = "outroot";
   Int_t whichPair = EncodePairType(kPip,kPim);
   if(argc>1) inDir = TString(argv[1]);
   if(argc>2) whichPair = (Int_t)strtof(argv[2],NULL);

   EventTree * ev = new EventTree(TString(inDir+"/*.root"),whichPair);

   TFile * outfile = new TFile("match.root","RECREATE");



   // define histograms
   const Int_t NBINS = 100;
  
     
   Float_t modulation,weight;

   printf("begin loop through %lld events...\n",ev->ENT);
   for(int i=0; i<ev->ENT; i++) {

     ev->GetEvent(i);

     if(ev->Valid()) {

       PhiHDist->Fill(ev->PhiH);
       PhiRDist->Fill(ev->PhiR);
       PhiHvsPhiR->Fill(ev->PhiR,ev->PhiH);

       PhiHR = Tools::AdjAngle(ev->PhiH - ev->PhiR);
       PhiHRDist->Fill(PhiHR);

       deltaPhi = Tools::AdjAngle(ev->hadPhi[qA] - ev->hadPhi[qB]);
       deltaPhiDist->Fill(deltaPhi);

       // set Asymmetry branches and fill modulation plots
       for(int m=0; m<Asymmetry::nMod; m++) {

         A[m]->Mh = ev->Mh;
         A[m]->x = ev->x;
         A[m]->z = ev->Zpair;
         A[m]->PhPerp = ev->PhPerp;
         A[m]->PhiH = ev->PhiH;
         A[m]->PhiR = ev->PhiR;
         //A[m]->PhiR = ev->PhiRq;

         A[m]->spinn = ev->SpinState();

         modulation = A[m]->EvalModulation();
         weight = A[m]->EvalWeight();

         ModVsZ[m]->Fill( ev->Zpair, modulation, weight );
         ModVsX[m]->Fill( ev->x, modulation, weight );
         ModVsMh[m]->Fill( ev->Mh, modulation, weight );
         ModVsPhPerp[m]->Fill( ev->PhPerp, modulation, weight );
         ModVsPhPerpMh[m]->Fill( ev->PhPerp/ev->Mh, modulation, weight );
         ModVsPhiH[m]->Fill( ev->PhiH, modulation, weight );
         ModVsPhiR[m]->Fill( ev->PhiR, modulation, weight );
       };

     };
   };

   deltaPhiDist->Write();
   PhiHDist->Write();
   PhiRDist->Write();
   PhiHvsPhiR->Write();
   PhiHRDist->Write();

   for(int m=0; m<Asymmetry::nMod; m++) {
     WriteCanvas(ModVsZ[m]);
     WriteCanvas(ModVsX[m]);
     WriteCanvas(ModVsMh[m]);
     WriteCanvas(ModVsPhPerp[m]);
     WriteCanvas(ModVsPhPerpMh[m]);
     WriteCanvas(ModVsPhiH[m]);
     WriteCanvas(ModVsPhiR[m]);
   };

   outfile->Close();

};



void WriteCanvas(TH2D * d) {
  TString canvName = Form("%s_canv",d->GetName());
  TCanvas * c = new TCanvas(canvName,canvName,1000,500);
  c->Divide(2,1);

  Float_t max = d->GetYaxis()->GetXmax();

  TProfile * dp = d->ProfileX();
  dp->SetLineColor(kBlack);
  dp->SetLineWidth(2);

  c->cd(1);
  c->GetPad(1)->SetLogz();
  d->Draw("colz");

  c->cd(2);
  dp->GetYaxis()->SetRangeUser(-max,max);
  dp->Draw();

  c->Write();
};

