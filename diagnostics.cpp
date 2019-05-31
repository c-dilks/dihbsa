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
TString corrTitle(TString var);
TString distTitle(TString var);

TString inDir;
Int_t whichPair;

int main(int argc, char** argv) {

   // ARGUMENTS
   inDir = "outroot";
   whichPair = pairPM;
   if(argc>1) inDir = TString(argv[1]);
   if(argc>2) whichPair = (Int_t)strtof(argv[2],NULL);

   EventTree * ev = new EventTree(TString(inDir+"/*.root"),whichPair);


   TFile * outfile = new TFile("plots.root","RECREATE");

   const Int_t NBINS = 100;
   Float_t deltaPhi;
   Float_t PhiHR;


   TH1F * WDist = new TH1F("WDist","W distribution (w/o W cut);W",NBINS,0,6);
   TH1F * XDist = new TH1F("XDist","x distribution;x",NBINS,0,1);
   TH2F * Q2vsW = new TH2F("Q2vsW","Q^{2} vs. W (w/o W cut);W;Q^{2}",
                                   NBINS,0,6,NBINS,0,12);
   TH2F * Q2vsX = new TH2F("Q2vsX","Q^{2} vs. x;x;Q^{2}",NBINS,0,1,NBINS,0,12);
   TH1F * YDist = new TH1F("YDist","y distribution (w/o y cut);y",NBINS,0,1);

   TH2F * hadECorr = new TH2F("hadECorr",corrTitle("E"),NBINS,0,10,NBINS,0,10);
   TH2F * hadPCorr = new TH2F("hadPCorr",corrTitle("p"),NBINS,0,10,NBINS,0,10);
   TH2F * hadPtCorr = new TH2F("hadPtCorr",corrTitle("p_{T}"),NBINS,0,4,NBINS,0,4);
   TH2F * hadEtaCorr = new TH2F("hadEtaCorr",corrTitle("#eta"),NBINS,0,5,NBINS,0,5);
   TH2F * hadPhiCorr = new TH2F("hadPhiCorr",corrTitle("#phi"),
                                             NBINS,-PI-1,PI+1,NBINS,-PI-1,PI+1);
   TH2F * hadZCorr = new TH2F("hadZCorr",corrTitle("z"),NBINS,0,1,NBINS,0,1);
   
   TH1F * hadEDist[2];
   TH1F * hadPDist[2];
   TH1F * hadPtDist[2];
   TH1F * hadEtaDist[2];
   TH1F * hadPhiDist[2];
   TH1F * hadZDist[2];
   for(int h=0; h<2; h++) {
     hadEDist[h] = new TH1F(TString(pmName(whichPair,h)+"hadEDist"),distTitle("E"),
       NBINS,0,10);
     hadPDist[h] = new TH1F(TString(pmName(whichPair,h)+"hadPDist"),distTitle("p"),
       NBINS,0,10);
     hadPtDist[h] = new TH1F(TString(pmName(whichPair,h)+"hadPtDist"),distTitle("p_{T}"),
       NBINS,0,4);
     hadEtaDist[h] = new TH1F(TString(pmName(whichPair,h)+"hadEtaDist"),distTitle("#eta"),
       NBINS,0,5);
     hadPhiDist[h] = new TH1F(TString(pmName(whichPair,h)+"hadPhiDist"),distTitle("#phi"),
       NBINS,-PI-1,PI+1);
     hadZDist[h] = new TH1F(TString(pmName(whichPair,h)+"hadZDist"),distTitle("z"),
       NBINS,0,1);
   };


   TString plotTitle = "#Delta#phi = #phi(" + pmTitle(whichPair,hP) + ")" +
                                 " - #phi(" + pmTitle(whichPair,hM) + 
                                 ") distribution;#Delta#phi";
   TH1F * deltaPhiDist = new TH1F("deltaPhiDist",plotTitle,NBINS,-PI-1,PI+1);

   TH1F * MhDist = new TH1F("MhDist","M_{h} distribution;M_{h}",NBINS,0,4);
   TH1F * ZpairDist = new TH1F("ZpairDist","z_{pair} distribution;z_{pair}",NBINS,0,1);
   TH1F * xFDist = new TH1F("xFDist","x_{F} distribution;x_{F}",NBINS,-2,2);
   TH1F * MmissDist = new TH1F("MmissDist","M_{miss} distribution;M_{miss}",NBINS,-2,6);
   
   TH1F * PhiHDist = new TH1F("PhiHDist","#phi_{h} distribution;#phi_{h}",
     NBINS,-PI-1,PI+1);
   TH1F * PhiRDist = new TH1F("PhiRDist","#phi_{R} distribution;#phi_{R}",
     NBINS,-PI-1,PI+1);
   TH2F * PhiHvsPhiR = new TH2F("PhiHvsPhiR","#phi_{h} vs. #phi_{R};#phi_{R};#phi_{h}",
     NBINS,-PI-1,PI+1,NBINS,-PI-1,PI+1);
   TH1F * PhiHRDist = new TH1F("PhiHRDist",
     "#phi_{h}-#phi_{R} distribution;#phi_{h}-#phi_{R}",
     NBINS,-PI-1,PI+1);

   plotTitle = "P_{h}^{perp}/M_{h} vs. sin(#phi_{h}-#phi_{R});";
   plotTitle += "sin(#phi_{h}-#phi_{R});P_{h}^{perp}/M_{h}";
   TH2F * g1perpWeightVsMod = new TH2F("g1perpWeightVsMod",plotTitle,
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



// make title for hadron correlation plots
TString corrTitle(TString var) {
  TString varX = var + "(" + pmTitle(whichPair,hM) + ")";
  TString varY = var + "(" + pmTitle(whichPair,hP) + ")";
  TString ret = varY + " vs. " + varX + ";" + varX + ";" + varY;
  return ret;
};


TString distTitle(TString var) {
  TString col[2]; 
  for(int cc=0; cc<2; cc++) {
    col[cc] = PartColorName(pmIdx(whichPair,cc)) + ":" + PartTitle(pmIdx(whichPair,cc));
  };
  TString ret = var + " distribution (" + col[0] + " " + col[1] + ");" + var;
  return ret;
};



// make canvas for hadron correlation plots
void HadronCompareCanv(TCanvas * canv, TH1F * dist[2], TH2F * corr) {

  for(int h=0; h<2; h++) {
    dist[h]->SetLineColor(PartColor(pmIdx(whichPair,h)));
    dist[h]->SetLineWidth(2);
  };

  canv->Divide(2,1);

  canv->cd(1);
  Int_t f = dist[hP]->GetMaximum() > dist[hM]->GetMaximum() ? hP:hM;
  dist[f]->Draw();
  dist[(f+1)%2]->Draw("SAME");

  canv->cd(2);
  corr->Draw("colz");
  canv->GetPad(2)->SetGrid(1,1);
};

