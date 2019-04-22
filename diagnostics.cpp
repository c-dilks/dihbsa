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

// DihBsa
#include "Constants.h"
#include "DIS.h"
#include "Trajectory.h"
#include "Dihadron.h"


int main(int argc, char** argv) {

  TChain * chain = new TChain("tree");
  chain->Add("outroot/out*.root");
  
  Int_t ENT = chain->GetEntries();
  //ENT = 10000;
  

   Float_t W,Q2,Nu,X,y;

   Float_t hadE[2]; // [enum plus_minus (0=+, 1=-)]
   Float_t hadP[2];
   Float_t hadPt[2];
   Float_t hadEta[2];
   Float_t hadPhi[2];

   Float_t Mh,Zpair,PhiH,PhiR,Mmiss,xF;
   Float_t Z[2];

   Int_t evnum,runnum;
   Int_t helicity;
   Float_t torus;
   Long64_t triggerBits;

   chain->SetBranchAddress("W",&W);
   chain->SetBranchAddress("Q2",&Q2);
   chain->SetBranchAddress("Nu",&Nu);
   chain->SetBranchAddress("X",&X);
   chain->SetBranchAddress("y",&y);

   chain->SetBranchAddress("hadE",hadE);
   chain->SetBranchAddress("hadP",hadP);
   chain->SetBranchAddress("hadPt",hadPt);
   chain->SetBranchAddress("hadEta",hadEta);
   chain->SetBranchAddress("hadPhi",hadPhi);

   chain->SetBranchAddress("Mh",&Mh);
   chain->SetBranchAddress("Z",Z);
   chain->SetBranchAddress("Zpair",&Zpair);
   chain->SetBranchAddress("PhiH",&PhiH);
   chain->SetBranchAddress("PhiR",&PhiR);
   chain->SetBranchAddress("Mmiss",&Mmiss);
   chain->SetBranchAddress("xF",&xF);

   chain->SetBranchAddress("runnum",&runnum);
   chain->SetBranchAddress("evnum",&evnum);
   chain->SetBranchAddress("helicity",&helicity);
   chain->SetBranchAddress("torus",&torus);
   chain->SetBranchAddress("triggerBits",&triggerBits);


   Bool_t cutQ2,cutW,cutY;
   Bool_t cutDihadron;

   Float_t deltaPhi;
   Float_t PhiHR;


   TFile * outfile = new TFile("plots.root","RECREATE");

   const Int_t NBINS = 100;

   TH1F * WDist = new TH1F("WDist","W distribution (w/o W cut);W",
     NBINS,0,6);
   TH2F * Q2vsW = new TH2F("Q2vsW","Q^{2} vs. W (w/o W cut);W;Q^{2}",
     NBINS,0,6,NBINS,0,10);
   TH2F * Q2vsX = new TH2F("Q2vsX","Q^{2} vs. x;x;Q^{2}",
     NBINS,0,1,NBINS,0,10);
   TH1F * YDist = new TH1F("YDist","y distribution (w/o y cut)",
     NBINS,0,1);

   TH2F * hadECorr = new TH2F("hadECorr",
     "E^{+} vs. E^{-};E^{-};E^{+}",
     NBINS,0,10,NBINS,0,10);
   TH2F * hadPCorr = new TH2F("hadPCorr",
     "p^{+} vs. p^{-};p^{-};p^{+}",
     NBINS,0,10,NBINS,0,10);
   TH2F * hadPtCorr = new TH2F("hadPtCorr",
     "p_{T}^{+} vs. p_{T}^{-};p_{T}^{-};p_{T}^{+}",
     NBINS,0,4,NBINS,0,4);
   TH2F * hadEtaCorr = new TH2F("hadEtaCorr",
     "#eta^{+} vs. #eta^{-};#eta^{-};#eta^{+}",
     NBINS,0,5,NBINS,0,5);
   TH2F * hadPhiCorr = new TH2F("hadPhiCorr",
     "#phi^{+} vs. #phi^{-};#phi^{-};#phi^{+}",
     NBINS,-PI-1,PI+1,NBINS,-PI-1,PI+1);
   TH2F * hadZCorr = new TH2F("hadZCorr",
     "z^{+} vs. z^{-};z^{-};z^{+}",
     NBINS,0,1,NBINS,0,1);
   
   TH1F * hadEDist[2];
   TH1F * hadPDist[2];
   TH1F * hadPtDist[2];
   TH1F * hadEtaDist[2];
   TH1F * hadPhiDist[2];
   TH1F * hadZDist[2];
   TString sgnSym[2];
   sgnSym[hP] = "+";
   sgnSym[hM] = "-";
   TString sgnTxt[2];
   sgnTxt[hP] = "p";
   sgnTxt[hM] = "m";
   TString plotTitle,plotName;
   for(int h=0; h<2; h++) {
     plotTitle = "E distribution (blue:#pi^{+} red:#pi^{-})";
     plotName = sgnTxt[h]+"hadEDist";
     hadEDist[h] = new TH1F(plotName,plotTitle,
       NBINS,0,10);
     plotTitle = "p distribution (blue:#pi^{+} red:#pi^{-})";
     plotName = sgnTxt[h]+"hadPDist";
     hadPDist[h] = new TH1F(plotName,plotTitle,
       NBINS,0,10);
     plotTitle = "p_{T} distribution (blue:#pi^{+} red:#pi^{-})";
     plotName = sgnTxt[h]+"hadPtDist";
     hadPtDist[h] = new TH1F(plotName,plotTitle,
       NBINS,0,4);
     plotTitle = "#eta distribution (blue:#pi^{+} red:#pi^{-})";
     plotName = sgnTxt[h]+"hadEtaDist";
     hadEtaDist[h] = new TH1F(plotName,plotTitle,
       NBINS,0,5);
     plotTitle = "#phi distribution (blue:#pi^{+} red:#pi^{-})";
     plotName = sgnTxt[h]+"hadPhiDist";
     hadPhiDist[h] = new TH1F(plotName,plotTitle,
       NBINS,-PI-1,PI+1);
     plotTitle = "Z distribution (blue:#pi^{+} red:#pi^{-})";
     plotName = sgnTxt[h]+"hadZDist";
     hadZDist[h] = new TH1F(plotName,plotTitle,
       NBINS,0,1);
   };
   Color_t plotColor[2];
   plotColor[hP] = kBlue;
   plotColor[hM] = kRed;
   for(int h=0; h<2; h++) {
     hadEDist[h]->SetLineColor(plotColor[h]);
     hadPDist[h]->SetLineColor(plotColor[h]);
     hadPtDist[h]->SetLineColor(plotColor[h]);
     hadEtaDist[h]->SetLineColor(plotColor[h]);
     hadPhiDist[h]->SetLineColor(plotColor[h]);
     hadZDist[h]->SetLineColor(plotColor[h]);

     hadEDist[h]->SetLineWidth(2);
     hadPDist[h]->SetLineWidth(2);
     hadPtDist[h]->SetLineWidth(2);
     hadEtaDist[h]->SetLineWidth(2);
     hadPhiDist[h]->SetLineWidth(2);
     hadZDist[h]->SetLineWidth(2);
   };



   TH1F * deltaPhiDist = new TH1F("deltaPhiDist",
     "#Delta#phi=#phi^{+}-#phi^{-} distribution;#Delta#phi",
     NBINS,-PI-1,PI+1);

   TH1F * MhDist = new TH1F("MhDist",
     "M_{h} distribution",
     NBINS,0,4);
   TH1F * ZpairDist = new TH1F("ZpairDist",
     "Z_{pair} distribution;Z_{pair}",
     NBINS,0,1);
   TH1F * xFDist = new TH1F("xFDist",
     "x_{F} distribution;x_{F}",
     NBINS,-2,2);
   TH1F * MmissDist = new TH1F("MmissDist",
     "M_{miss} distribution;M_{miss}",
     NBINS,-2,6);
   
   TH1F * PhiHDist = new TH1F("PhiHDist",
     "#phi_{h} distribution;#phi_{h}",
     NBINS,-PI-1,PI+1);
   TH1F * PhiRDist = new TH1F("PhiRDist",
     "#phi_{R} distribution;#phi_{R}",
     NBINS,-PI-1,PI+1);
   TH2F * PhiHvsPhiR = new TH2F("PhiHvsPhiR",
     "#phi_{h} vs. #phi_{R};#phi_{R};#phi_{h}",
     NBINS,-PI-1,PI+1,
     NBINS,-PI-1,PI+1);
   TH1F * PhiHRDist = new TH1F("PhiHRDist",
     "#phi_{h}-#phi_{R} distribution;#phi_{h}-#phi_{R}",
     NBINS,-PI-1,PI+1);








   printf("begin loop through %d events...\n",ENT);
   for(int i=0; i<ENT; i++) {
     chain->GetEntry(i);

     if(i%10000==0) printf("%.2f%%\n",100*(float)i/((float)ENT));


     cutQ2 = Q2 > 1.0;
     cutW = W > 2.0;
     cutY = y < 0.8;


     cutDihadron = true;
     cutDihadron = cutDihadron && Z[hP] > 0.1 && Z[hM] > 0.1;
     cutDihadron = cutDihadron && Zpair < 0.95;
     cutDihadron = cutDihadron && Mmiss > 1.05;
     cutDihadron = cutDihadron && xF > 0;
     cutDihadron = cutDihadron && hadP[hP] > 1.0 && hadP[hM] > 1.0;

     if(!cutDihadron) continue;


     if(cutQ2 && cutY) {
       WDist->Fill(W);
       Q2vsW->Fill(W,Q2);
     };

     if(cutQ2 && cutW) {
       YDist->Fill(y);
     };


     if(cutQ2 && cutW && cutY) {
       Q2vsX->Fill(X,Q2);


       hadECorr->Fill(hadE[hM],hadE[hP]);
       hadPCorr->Fill(hadP[hM],hadP[hP]);
       hadPtCorr->Fill(hadPt[hM],hadPt[hP]);
       hadEtaCorr->Fill(hadEta[hM],hadEta[hP]);
       hadPhiCorr->Fill(hadPhi[hM],hadPhi[hP]);
       hadZCorr->Fill(Z[hM],Z[hP]);

       for(int h=0; h<2; h++) {
         hadEDist[h]->Fill(hadE[h]);
         hadPDist[h]->Fill(hadP[h]);
         hadPtDist[h]->Fill(hadPt[h]);
         hadEtaDist[h]->Fill(hadEta[h]);
         hadPhiDist[h]->Fill(hadPhi[h]);
         hadZDist[h]->Fill(Z[h]);
       };

       deltaPhi = hadPhi[hP] - hadPhi[hM];
       while(deltaPhi>PI) deltaPhi-=2*PI;
       while(deltaPhi<-PI) deltaPhi+=2*PI;
       deltaPhiDist->Fill(deltaPhi);

       MhDist->Fill(Mh);
       ZpairDist->Fill(Zpair);
       xFDist->Fill(xF);
       MmissDist->Fill(Mmiss);

       PhiHDist->Fill(PhiH);
       PhiRDist->Fill(PhiR);
       PhiHvsPhiR->Fill(PhiR,PhiH);

       PhiHR = PhiH-PhiR;
       while(PhiHR>PI) PhiHR-=2*PI;
       while(PhiHR<-PI) PhiHR+=2*PI;
       PhiHRDist->Fill(PhiHR);
     };


   };



   Int_t f;
   TCanvas * hadECanv = new TCanvas("hadECanv","hadECanv",
     1000,800);
   hadECanv->Divide(2,1);
   hadECanv->cd(1);
   f = hadEDist[hP]->GetMaximum() > hadEDist[hM]->GetMaximum() ?
     hP:hM;
   hadEDist[f]->Draw();
   hadEDist[(f+1)%2]->Draw("SAME");
   hadECanv->cd(2);
   hadECorr->Draw("colz");
   TCanvas * hadPCanv = new TCanvas("hadPCanv","hadPCanv",
     1000,800);
   hadPCanv->Divide(2,1);
   hadPCanv->cd(1);
   f = hadPDist[hP]->GetMaximum() > hadPDist[hM]->GetMaximum() ?
     hP:hM;
   hadPDist[f]->Draw();
   hadPDist[(f+1)%2]->Draw("SAME");
   hadPCanv->cd(2);
   hadPCorr->Draw("colz");
   TCanvas * hadPtCanv = new TCanvas("hadPtCanv","hadPtCanv",
     1000,800);
   hadPtCanv->Divide(2,1);
   hadPtCanv->cd(1);
   f = hadPtDist[hP]->GetMaximum() > hadPtDist[hM]->GetMaximum() ?
     hP:hM;
   hadPtDist[f]->Draw();
   hadPtDist[(f+1)%2]->Draw("SAME");
   hadPtCanv->cd(2);
   hadPtCorr->Draw("colz");
   TCanvas * hadEtaCanv = new TCanvas("hadEtaCanv","hadEtaCanv",
     1000,800);
   hadEtaCanv->Divide(2,1);
   hadEtaCanv->cd(1);
   f = hadEtaDist[hP]->GetMaximum() > hadEtaDist[hM]->GetMaximum() ?
     hP:hM;
   hadEtaDist[f]->Draw();
   hadEtaDist[(f+1)%2]->Draw("SAME");
   hadEtaCanv->cd(2);
   hadEtaCorr->Draw("colz");
   TCanvas * hadPhiCanv = new TCanvas("hadPhiCanv","hadPhiCanv",
     1000,800);
   hadPhiCanv->Divide(2,1);
   hadPhiCanv->cd(1);
   f = hadPhiDist[hP]->GetMaximum() > hadPhiDist[hM]->GetMaximum() ?
     hP:hM;
   hadPhiDist[f]->Draw();
   hadPhiDist[(f+1)%2]->Draw("SAME");
   hadPhiCanv->cd(2);
   hadPhiCorr->Draw("colz");
   TCanvas * hadZCanv = new TCanvas("hadZCanv","hadZCanv",
     1000,800);
   hadZCanv->Divide(2,1);
   hadZCanv->cd(1);
   f = hadZDist[hP]->GetMaximum() > hadZDist[hM]->GetMaximum() ?
     hP:hM;
   hadZDist[f]->Draw();
   hadZDist[(f+1)%2]->Draw("SAME");
   hadZCanv->cd(2);
   hadZCorr->Draw("colz");


   WDist->Write();
   Q2vsW->Write();
   Q2vsX->Write();
   YDist->Write();

   hadECanv->Write();
   hadPCanv->Write();
   hadPtCanv->Write();
   hadEtaCanv->Write();
   hadPhiCanv->Write();
   hadZCanv->Write();

   deltaPhiDist->Write();

   MhDist->Write();
   ZpairDist->Write();
   xFDist->Write();
   MmissDist->Write();

   PhiHDist->Write();
   PhiRDist->Write();
   PhiHvsPhiR->Write();

   PhiHRDist->Write();

};
  

