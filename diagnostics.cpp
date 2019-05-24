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
#include "Tools.h"
#include "DIS.h"
#include "Trajectory.h"
#include "Dihadron.h"
#include "EventTree.h"

void HadronCompareCanv(TCanvas * canv, TH1F * dist[2], TH2F * corr);

int main(int argc, char** argv) {

   // ARGUMENTS
   TString inDir = "outroot";
   Int_t pairSetting = pairPM;
   if(argc>1) inDir = TString(argv[1]);
   if(argc>2) pairSetting = (Int_t)strtof(argv[2],NULL);

   EventTree * ev = new EventTree(TString(inDir+"/*.root"));


   TFile * outfile = new TFile("plots.root","RECREATE");

   const Int_t NBINS = 100;
   Float_t deltaPhi;
   Float_t PhiHR;

   TString plotTitle,plotName;

   TH1F * WDist = new TH1F("WDist","W distribution (w/o W cut);W",
     NBINS,0,6);
   TH1F * XDist = new TH1F("XDist","x distribution",
     NBINS,0,1);
   TH2F * Q2vsW = new TH2F("Q2vsW","Q^{2} vs. W (w/o W cut);W;Q^{2}",
     NBINS,0,6,NBINS,0,12);
   TH2F * Q2vsX = new TH2F("Q2vsX","Q^{2} vs. x;x;Q^{2}",
     NBINS,0,1,NBINS,0,12);
   TH1F * YDist = new TH1F("YDist","y distribution (w/o y cut)",
     NBINS,0,1);

   plotTitle = Form("E^{%s} vs. E^{%s};E^{%s};E^{%s}",
     PMsym(pairSetting,hP).Data(),PMsym(pairSetting,hM),
     PMsym(pairSetting,hM).Data(),PMsym(pairSetting,hP));
   TH2F * hadECorr = new TH2F("hadECorr",
     plotTitle.Data(),
     NBINS,0,10,NBINS,0,10);
   plotTitle = Form("p^{%s} vs. p^{%s};p^{%s};p^{%s}",
     PMsym(pairSetting,hP).Data(),PMsym(pairSetting,hM),
     PMsym(pairSetting,hM).Data(),PMsym(pairSetting,hP));
   TH2F * hadPCorr = new TH2F("hadPCorr",
     plotTitle.Data(),
     NBINS,0,10,NBINS,0,10);
   plotTitle = Form("p_{T}^{%s} vs. p_{T}^{%s};p_{T}^{%s};p_{T}^{%s}",
     PMsym(pairSetting,hP).Data(),PMsym(pairSetting,hM),
     PMsym(pairSetting,hM).Data(),PMsym(pairSetting,hP));
   TH2F * hadPtCorr = new TH2F("hadPtCorr",
     plotTitle.Data(),
     NBINS,0,4,NBINS,0,4);
   plotTitle = Form("#eta^{%s} vs. #eta^{%s};#eta^{%s};#eta^{%s}",
     PMsym(pairSetting,hP).Data(),PMsym(pairSetting,hM),
     PMsym(pairSetting,hM).Data(),PMsym(pairSetting,hP));
   TH2F * hadEtaCorr = new TH2F("hadEtaCorr",
     plotTitle.Data(),
     NBINS,0,5,NBINS,0,5);
   plotTitle = Form("#phi^{%s} vs. #phi^{%s};#phi^{%s};#phi^{%s}",
     PMsym(pairSetting,hP).Data(),PMsym(pairSetting,hM),
     PMsym(pairSetting,hM).Data(),PMsym(pairSetting,hP));
   TH2F * hadPhiCorr = new TH2F("hadPhiCorr",
     plotTitle.Data(),
     NBINS,-PI-1,PI+1,NBINS,-PI-1,PI+1);
   plotTitle = Form("z^{%s} vs. z^{%s};z^{%s};z^{%s}",
     PMsym(pairSetting,hP).Data(),PMsym(pairSetting,hM),
     PMsym(pairSetting,hM).Data(),PMsym(pairSetting,hP));
   TH2F * hadZCorr = new TH2F("hadZCorr",
     plotTitle.Data(),
     NBINS,0,1,NBINS,0,1);
   
   TH1F * hadEDist[2];
   TH1F * hadPDist[2];
   TH1F * hadPtDist[2];
   TH1F * hadEtaDist[2];
   TH1F * hadPhiDist[2];
   TH1F * hadZDist[2];
   for(int h=0; h<2; h++) {
     plotTitle = Form("E distribution (blue:#pi^{%s} red:#pi^{%s})",
       PMsym(pairSetting,hP).Data(),PMsym(pairSetting,hM).Data());
     plotName = PMstr(pairSetting,h)+"hadEDist";
     hadEDist[h] = new TH1F(plotName,plotTitle,
       NBINS,0,10);
     plotTitle = Form("p distribution (blue:#pi^{%s} red:#pi^{%s})",
       PMsym(pairSetting,hP).Data(),PMsym(pairSetting,hM).Data());
     plotName = PMstr(pairSetting,h)+"hadPDist";
     hadPDist[h] = new TH1F(plotName,plotTitle,
       NBINS,0,10);
     plotTitle = Form("p_{T} distribution (blue:#pi^{%s} red:#pi^{%s})",
       PMsym(pairSetting,hP).Data(),PMsym(pairSetting,hM).Data());
     plotName = PMstr(pairSetting,h)+"hadPtDist";
     hadPtDist[h] = new TH1F(plotName,plotTitle,
       NBINS,0,4);
     plotTitle = Form("#eta distribution (blue:#pi^{%s} red:#pi^{%s})",
       PMsym(pairSetting,hP).Data(),PMsym(pairSetting,hM).Data());
     plotName = PMstr(pairSetting,h)+"hadEtaDist";
     hadEtaDist[h] = new TH1F(plotName,plotTitle,
       NBINS,0,5);
     plotTitle = Form("#phi distribution (blue:#pi^{%s} red:#pi^{%s})",
       PMsym(pairSetting,hP).Data(),PMsym(pairSetting,hM).Data());
     plotName = PMstr(pairSetting,h)+"hadPhiDist";
     hadPhiDist[h] = new TH1F(plotName,plotTitle,
       NBINS,-PI-1,PI+1);
     plotTitle = Form("Z distribution (blue:#pi^{%s} red:#pi^{%s})",
       PMsym(pairSetting,hP).Data(),PMsym(pairSetting,hM).Data());
     plotName = PMstr(pairSetting,h)+"hadZDist";
     hadZDist[h] = new TH1F(plotName,plotTitle,
       NBINS,0,1);
   };



   plotTitle = Form("#Delta#phi=#phi^{%s}-#phi^{%s} distribution;#Delta#phi",
     PMsym(pairSetting,hP).Data(),PMsym(pairSetting,hM).Data());
   TH1F * deltaPhiDist = new TH1F("deltaPhiDist",
     plotTitle,
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

   TH2F * g1perpWeightVsMod = new TH2F("g1perpWeightVsMod",
     "P_{h}^{perp}/M_{h} vs. sin(#phi_{h}-#phi_{R});sin(#phi_{h}-#phi_{R});P_{h}^{perp}/M_{h}",
     NBINS,-1.1,1.1,
     NBINS,0,6);





   printf("begin loop through %lld events...\n",ev->ENT);
   for(int i=0; i<ev->ENT; i++) {

     ev->GetEvent(i);

     if(ev->cutDihadron) {


       if(ev->cutQ2 && ev->cutY) {
         WDist->Fill(ev->W);
         Q2vsW->Fill(ev->W,ev->Q2);
       };

       if(ev->cutQ2 && ev->cutW) {
         YDist->Fill(ev->y);
       };


       if(ev->cutQ2 && ev->cutW && ev->cutY) {

         Q2vsX->Fill(ev->x,ev->Q2);
         XDist->Fill(ev->x);


         hadECorr->Fill(ev->hadE[hM],ev->hadE[hP]);
         hadPCorr->Fill(ev->hadP[hM],ev->hadP[hP]);
         hadPtCorr->Fill(ev->hadPt[hM],ev->hadPt[hP]);
         hadEtaCorr->Fill(ev->hadEta[hM],ev->hadEta[hP]);
         hadPhiCorr->Fill(ev->hadPhi[hM],ev->hadPhi[hP]);
         hadZCorr->Fill(ev->Z[hM],ev->Z[hP]);

         for(int h=0; h<2; h++) {
           hadEDist[h]->Fill(ev->hadE[h]);
           hadPDist[h]->Fill(ev->hadP[h]);
           hadPtDist[h]->Fill(ev->hadPt[h]);
           hadEtaDist[h]->Fill(ev->hadEta[h]);
           hadPhiDist[h]->Fill(ev->hadPhi[h]);
           hadZDist[h]->Fill(ev->Z[h]);
         };

         deltaPhi = Tools::AdjAngle(ev->hadPhi[hP] - ev->hadPhi[hM]);
         deltaPhiDist->Fill(deltaPhi);

         MhDist->Fill(ev->Mh);
         ZpairDist->Fill(ev->Zpair);
         xFDist->Fill(ev->xF);
         MmissDist->Fill(ev->Mmiss);

         PhiHDist->Fill(ev->PhiH);
         PhiRDist->Fill(ev->PhiR);
         PhiHvsPhiR->Fill(ev->PhiR,ev->PhiH);

         PhiHR = Tools::AdjAngle(ev->PhiH - ev->PhiR);
         PhiHRDist->Fill(PhiHR);
         g1perpWeightVsMod->Fill(TMath::Sin(PhiHR),ev->PhPerp/ev->Mh);
       };


     };
   };


   WDist->Write();
   XDist->Write();
   Q2vsW->Write();
   Q2vsX->Write();
   YDist->Write();

   TCanvas * hadECanv = new TCanvas("hadECanv","hadECanv",1000,800);
   TCanvas * hadPCanv = new TCanvas("hadPCanv","hadPCanv",1000,800);
   TCanvas * hadPtCanv = new TCanvas("hadPtCanv","hadPtCanv",1000,800);
   TCanvas * hadEtaCanv = new TCanvas("hadEtaCanv","hadEtaCanv",1000,800);
   TCanvas * hadPhiCanv = new TCanvas("hadPhiCanv","hadPhiCanv",1000,800);
   TCanvas * hadZCanv = new TCanvas("hadZCanv","hadZCanv",1000,800);

   HadronCompareCanv(hadECanv, hadEDist, hadECorr);
   HadronCompareCanv(hadPCanv, hadPDist, hadPCorr);
   HadronCompareCanv(hadPtCanv, hadPtDist, hadPtCorr);
   HadronCompareCanv(hadEtaCanv, hadEtaDist, hadEtaCorr);
   HadronCompareCanv(hadPhiCanv, hadPhiDist, hadPhiCorr);
   HadronCompareCanv(hadZCanv, hadZDist, hadZCorr);

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
   g1perpWeightVsMod->Write();

   outfile->Close();

};



void HadronCompareCanv(TCanvas * canv, TH1F * dist[2], TH2F * corr) {

  dist[hP]->SetLineColor(kBlue);
  dist[hM]->SetLineColor(kRed);
  for(int h=0; h<2; h++) dist[h]->SetLineWidth(2);

  canv->Divide(2,1);

  canv->cd(1);
  Int_t f = dist[hP]->GetMaximum() > dist[hM]->GetMaximum() ? hP:hM;
  dist[f]->Draw();
  dist[(f+1)%2]->Draw("SAME");

  canv->cd(2);
  corr->Draw("colz");
  canv->GetPad(2)->SetGrid(1,1);
};

