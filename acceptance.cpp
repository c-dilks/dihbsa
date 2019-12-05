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
#include "Trajectory.h"
#include "Dihadron.h"
#include "EventTree.h"
#include "Binning.h"
#include "Asymmetry.h"

void WriteCanvas(TH2D * d);

int main(int argc, char** argv) {

   // ARGUMENTS
   TString inDir = "outroot";
   Int_t whichPair = EncodePairType(kPip,kPim);
   if(argc>1) inDir = TString(argv[1]);
   if(argc>2) whichPair = (Int_t)strtof(argv[2],NULL);

   EventTree * ev = new EventTree(TString(inDir+"/*.root"),whichPair);


   TFile * outfile = new TFile("acc.root","RECREATE");

   const Int_t NBINS = 100;
   Float_t deltaPhi;
   Float_t PhiHR;


   // instantiate Asymmetry array, used for calculating azimuthal modulations
   Binning * B = new Binning(whichPair);
   Asymmetry * A[Asymmetry::nMod];
   for(int m=0; m<Asymmetry::nMod; m++) {
     B->AsymModulation = m;
     A[m] = new Asymmetry(B);
     printf("m=%d  %s\n",m,A[m]->ModulationTitle.Data());
   };

  
   TH1F * deltaPhiDist = new TH1F("deltaPhiDist",
     "#Delta#phi=#phi^{+}-#phi^{-} distribution;#Delta#phi",
     NBINS,-PI,PI);
   
   TH1F * PhiHDist = new TH1F("PhiHDist",
     "#phi_{h} distribution;#phi_{h}",
     NBINS,-PI,PI);
   TH1F * PhiRDist = new TH1F("PhiRDist",
     "#phi_{R} distribution;#phi_{R}",
     NBINS,-PI,PI);
   TH2F * PhiHvsPhiR = new TH2F("PhiHvsPhiR",
     "#phi_{h} vs. #phi_{R};#phi_{R};#phi_{h}",
     NBINS,-PI,PI,
     NBINS,-PI,PI);

   TH1F * PhiHRDist = new TH1F("PhiHRDist",
     "#phi_{h}-#phi_{R} distribution;#phi_{h}-#phi_{R}",
     NBINS,-PI,PI);


   TH2D * ModVsZ[Asymmetry::nMod];
   TH2D * ModVsX[Asymmetry::nMod];
   TH2D * ModVsMh[Asymmetry::nMod];
   TH2D * ModVsPhPerp[Asymmetry::nMod];
   TH2D * ModVsPhPerpMh[Asymmetry::nMod];
   TH2D * ModVsPhiH[Asymmetry::nMod];
   TH2D * ModVsPhiR[Asymmetry::nMod];
   TString plotTitle,plotName;

   for(int m=0; m<Asymmetry::nMod; m++) {
     ModVsZ[m] = new TH2D(
       TString(A[m]->ModulationName + "_vs_Z"),
       TString(plotTitle = A[m]->ModulationTitle + " vs. z"),
       NBINS,0,1,
       NBINS,-1*(A[m]->modMax),A[m]->modMax);
     ModVsX[m] = new TH2D(
       TString(A[m]->ModulationName + "_vs_X"),
       TString(plotTitle = A[m]->ModulationTitle + " vs. x"),
       NBINS,0,1,
       NBINS,-1*(A[m]->modMax),A[m]->modMax);
     ModVsMh[m] = new TH2D(
       TString(A[m]->ModulationName + "_vs_Mh"),
       TString(plotTitle = A[m]->ModulationTitle + " vs. M_{h}"),
       NBINS,0,3,
       NBINS,-1*(A[m]->modMax),A[m]->modMax);
     ModVsPhPerp[m] = new TH2D(
       TString(A[m]->ModulationName + "_vs_PhPerp"),
       TString(plotTitle = A[m]->ModulationTitle + " vs. P_{h}^{perp}"),
       NBINS,0,2,
       NBINS,-1*(A[m]->modMax),A[m]->modMax);
     ModVsPhPerpMh[m] = new TH2D(
       TString(A[m]->ModulationName + "_vs_PhPerpMh"),
       TString(plotTitle = A[m]->ModulationTitle + " vs. P_{h}^{perp}/M_{h}"),
       NBINS,0,5,
       NBINS,-1*(A[m]->modMax),A[m]->modMax);
     ModVsPhiH[m] = new TH2D(
       TString(A[m]->ModulationName + "_vs_PhiH"),
       TString(plotTitle = A[m]->ModulationTitle + " vs. #phi_{h}"),
       NBINS,-PI,PI,
       NBINS,-1*(A[m]->modMax),A[m]->modMax);
     ModVsPhiR[m] = new TH2D(
       TString(A[m]->ModulationName + "_vs_PhiR"),
       TString(plotTitle = A[m]->ModulationTitle + " vs. #phi_{R}"),
       NBINS,-PI,PI,
       NBINS,-1*(A[m]->modMax),A[m]->modMax);
   };

     
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

