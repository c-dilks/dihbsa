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

void HadronCompareCanv(TCanvas * canv, TH1D * dist[2], TH2D * corr);
TString corrTitle(TString var);
TString distTitle(TString var);

TString inDir;
Int_t whichPair;
Int_t whichHad[2];
TString hadName[2];
TString hadTitle[2];

int main(int argc, char** argv) {

   // ARGUMENTS
   inDir = "outroot";
   whichPair = Tools::EncodePairType(kPip,kPim);
   if(argc>1) inDir = TString(argv[1]);
   if(argc>2) whichPair = (Int_t)strtof(argv[2],NULL);
   
   // get hadron pair from whichPair; note that in the print out, the 
   // order of hadron 0 and 1 is set by Constants::dihHadIdx
   printf("whichPair = 0x%x\n",whichPair);
   Tools::DecodePairType(whichPair,whichHad[qA],whichHad[qB]);
   for(int h=0; h<2; h++) {
     hadName[h] = PairHadName(whichHad[qA],whichHad[qB],h);
     hadTitle[h] = PairHadTitle(whichHad[qA],whichHad[qB],h);
     printf("hadron %d:  idx=%d  name=%s  title=%s\n",
       h,dihHadIdx(whichHad[qA],whichHad[qB],h),hadName[h].Data(),hadTitle[h].Data());
   };

   EventTree * ev = new EventTree(TString(inDir+"/*.root"),whichPair);


   TFile * outfile = new TFile("plots.root","RECREATE");

   const Int_t NBINS = 100;
   Float_t deltaPhi;
   Float_t PhiHR;


   // DIS kinematics
   TH1D * WDist = new TH1D("WDist","W distribution (w/o W cut);W",NBINS,0,6);
   TH1D * XDist = new TH1D("XDist","x distribution;x",NBINS,0,1);
   TH2D * Q2vsW = new TH2D("Q2vsW","Q^{2} vs. W (w/o W cut);W;Q^{2}",
                                   NBINS,0,6,NBINS,0,12);
   TH2D * Q2vsX = new TH2D("Q2vsX","Q^{2} vs. x;x;Q^{2}",NBINS,0,1,NBINS,0,12);
   TH1D * YDist = new TH1D("YDist","y distribution (w/o y cut);y",NBINS,0,1);
   
   // electron kinematics
   TH1D * eleEDist = new TH1D("eleEDist","e^{-} E distribution",NBINS,0,12);
   TH1D * elePtDist = new TH1D("elePtDist","e^{-} p_{T} distribution",NBINS,0,4);
   TH1D * eleEtaDist = new TH1D("eleEtaDist","e^{-} #phi distribution",NBINS,-3,6);
   TH1D * elePhiDist = new TH1D("elePhiDist","e^{-} #phi distribution",NBINS,-PIe,PIe);
   TH2D * elePtVsPhi = new TH2D("elePtvsPhi","e^{-} p_{T} vs #phi;#phi;#p_{T}",
     NBINS,-PIe,PIe,NBINS,0,4);

   // dihadron's hadron kinematic correlations
   TH2D * hadECorr = new TH2D("hadECorr",corrTitle("E"),NBINS,0,10,NBINS,0,10);
   TH2D * hadPCorr = new TH2D("hadPCorr",corrTitle("p"),NBINS,0,10,NBINS,0,10);
   TH2D * hadPtCorr = new TH2D("hadPtCorr",corrTitle("p_{T}"),NBINS,0,4,NBINS,0,4);
   TH2D * hadEtaCorr = new TH2D("hadEtaCorr",corrTitle("#eta"),NBINS,0,5,NBINS,0,5);
   TH2D * hadPhiCorr = new TH2D("hadPhiCorr",corrTitle("#phi"),
                                             NBINS,-PIe,PIe,NBINS,-PIe,PIe);
   TH2D * hadZCorr = new TH2D("hadZCorr",corrTitle("z"),NBINS,0,1,NBINS,0,1);
   
   // dihadron's hadron kinematics
   TH1D * hadEDist[2];
   TH1D * hadPDist[2];
   TH1D * hadPtDist[2];
   TH1D * hadEtaDist[2];
   TH1D * hadPhiDist[2];
   TH1D * hadZDist[2];
   for(int h=0; h<2; h++) {
     hadEDist[h] = new TH1D(TString(hadName[h]+"hadEDist"),distTitle("E"),
       NBINS,0,10);
     hadPDist[h] = new TH1D(TString(hadName[h]+"hadPDist"),distTitle("p"),
       NBINS,0,10);
     hadPtDist[h] = new TH1D(TString(hadName[h]+"hadPtDist"),distTitle("p_{T}"),
       NBINS,0,4);
     hadEtaDist[h] = new TH1D(TString(hadName[h]+"hadEtaDist"),distTitle("#eta"),
       NBINS,0,5);
     hadPhiDist[h] = new TH1D(TString(hadName[h]+"hadPhiDist"),distTitle("#phi"),
       NBINS,-PIe,PIe);
     hadZDist[h] = new TH1D(TString(hadName[h]+"hadZDist"),distTitle("z"),
       NBINS,0,1);
   };


   // dihadron kinematics
   TString plotTitle = "#Delta#phi = #phi(" + hadTitle[qA] + ")" +
                                 " - #phi(" + hadTitle[qB] + 
                                 ") distribution;#Delta#phi";
   TH1D * deltaPhiDist = new TH1D("deltaPhiDist",plotTitle,NBINS,-PIe,PIe);

   TH1D * MhDist = new TH1D("MhDist","M_{h} distribution;M_{h}",NBINS,0,4);
   TH1D * ZpairDist = new TH1D("ZpairDist","z_{pair} distribution;z_{pair}",NBINS,0,1);
   TH1D * xFDist = new TH1D("xFDist","x_{F} distribution;x_{F}",NBINS,-2,2);
   TH1D * MmissDist = new TH1D("MmissDist","M_{X} distribution;M_{X}",NBINS,-2,6);
   
   TH1D * PhiHDist = new TH1D("PhiHDist","#phi_{h} distribution;#phi_{h}",
     NBINS,-PIe,PIe);
   TH1D * PhiRDist = new TH1D("PhiRDist","#phi_{R} distribution;#phi_{R}",
     NBINS,-PIe,PIe);
   TH2D * PhiHvsPhiR = new TH2D("PhiHvsPhiR","#phi_{h} vs. #phi_{R};#phi_{R};#phi_{h}",
     NBINS,-PIe,PIe,NBINS,-PIe,PIe);
   TH1D * PhiHRDist = new TH1D("PhiHRDist",
     "#phi_{h}-#phi_{R} distribution;#phi_{h}-#phi_{R}",
     NBINS,-PIe,PIe);

   plotTitle = "P_{h}^{perp}/M_{h} vs. sin(#phi_{h}-#phi_{R});";
   plotTitle += "sin(#phi_{h}-#phi_{R});P_{h}^{perp}/M_{h}";
   TH2D * g1perpWeightVsMod = new TH2D("g1perpWeightVsMod",plotTitle,
     NBINS,-1.1,1.1,
     NBINS,0,6);


   // hadron type matrix
   TH2D * hadTypeMatrix = new TH2D("hadTypeMatrix",
     "Dihadron hadron types matrix;hadron 2;hadron 1",
     nObservables,0,nObservables,nObservables,0,nObservables);
   for(int h1=0; h1<nObservables; h1++) {
     for(int h2=0; h2<nObservables; h2++) {
       hadTypeMatrix->GetYaxis()->SetBinLabel(h1+1,ObsTitle(h1));
       hadTypeMatrix->GetXaxis()->SetBinLabel(h2+1,ObsTitle(h2));
     };
   };


   // fiducial phi mask
   TH1D * fiducialPhiMask = new TH1D("fiducialPhiMask",
     "#phi fiducial regions",NBINS,-PIe,PIe);





   printf("begin loop through %lld events...\n",ev->ENT);
   Int_t hadI[2];
   for(int i=0; i<ev->ENT; i++) {

     ev->GetEvent(i);

     // fill hadron types matrix (note cutDihadronKinematics does all dihadron
     // kinematic cuts, except for demanding the hadron Idx's are the requested ones)
     for(int h=0; h<2; h++) hadI[h] = IO(ev->hadIdx[h]);
     if(ev->cutDihadronKinematics && ev->cutQ2 && ev->cutW && ev->cutY) {
       hadTypeMatrix->Fill(hadI[qB],hadI[qA]);
       // fill transpose elements too (makes symmetric matrix), but don't double-fill
       // diagonal elements
       if(hadI[qA]!=hadI[qB]) hadTypeMatrix->Fill(hadI[qA],hadI[qB]);
     };


     // fill DIS kinematic plots
     if(ev->cutDihadron) {

       if(ev->cutQ2 && ev->cutY) {
         WDist->Fill(ev->W);
         Q2vsW->Fill(ev->W,ev->Q2);
       };

       if(ev->cutQ2 && ev->cutW) {
         YDist->Fill(ev->y);
       };
     };


     // fill dihadron kinematics plots
     if(ev->cutDihadron && ev->cutQ2 && ev->cutW && ev->cutY) {

       eleEDist->Fill(ev->eleE);
       elePtDist->Fill(ev->elePt);
       eleEtaDist->Fill(ev->eleEta);
       elePhiDist->Fill(ev->elePhi);
       elePtVsPhi->Fill(ev->elePhi,ev->elePt);

       if(Tools::PhiFiducialCut(ev->elePhi)) fiducialPhiMask->Fill(ev->elePhi);


       Q2vsX->Fill(ev->x,ev->Q2);
       XDist->Fill(ev->x);


       hadECorr->Fill(ev->hadE[qB],ev->hadE[qA]);
       hadPCorr->Fill(ev->hadP[qB],ev->hadP[qA]);
       hadPtCorr->Fill(ev->hadPt[qB],ev->hadPt[qA]);
       hadEtaCorr->Fill(ev->hadEta[qB],ev->hadEta[qA]);
       hadPhiCorr->Fill(ev->hadPhi[qB],ev->hadPhi[qA]);
       hadZCorr->Fill(ev->Z[qB],ev->Z[qA]);

       for(int h=0; h<2; h++) {
         hadEDist[h]->Fill(ev->hadE[h]);
         hadPDist[h]->Fill(ev->hadP[h]);
         hadPtDist[h]->Fill(ev->hadPt[h]);
         hadEtaDist[h]->Fill(ev->hadEta[h]);
         hadPhiDist[h]->Fill(ev->hadPhi[h]);
         hadZDist[h]->Fill(ev->Z[h]);
       };

       deltaPhi = Tools::AdjAngle(ev->hadPhi[qA] - ev->hadPhi[qB]);
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


   }; // eo event loop


   WDist->Write();
   XDist->Write();
   Q2vsW->Write();
   Q2vsX->Write();
   YDist->Write();

   eleEDist->Write();
   elePtDist->Write();
   eleEtaDist->Write();
   elePhiDist->Write();
   elePtVsPhi->Write();
   fiducialPhiMask->Write();



   hadTypeMatrix->Write();

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
  TString varX = var + "(" + hadTitle[qB] + ")";
  TString varY = var + "(" + hadTitle[qA] + ")";
  TString ret = varY + " vs. " + varX + ";" + varX + ";" + varY;
  return ret;
};


TString distTitle(TString var) {
  TString col[2]; 
  if(whichHad[qA]!=whichHad[qB]) {
    for(int cc=0; cc<2; cc++) {
      col[cc] = PartColorName(dihHadIdx(whichHad[qA],whichHad[qB],cc)) + 
                ":" + hadTitle[cc];
    };
  } else {
    col[qA] = "solid:" + hadTitle[qA];
    col[qB] = "dashed:" + hadTitle[qB]; 
  };
  TString ret = var + " distribution (" + col[qA] + " " + col[qB] + ");" + var;
  return ret;
};



// make canvas for hadron correlation plots
void HadronCompareCanv(TCanvas * canv, TH1D * dist[2], TH2D * corr) {

  for(int h=0; h<2; h++) {
    dist[h]->SetLineColor(PartColor(dihHadIdx(whichHad[qA],whichHad[qB],h)));
    dist[h]->SetLineWidth(2);
  };
  if(whichHad[qA]==whichHad[qB]) dist[qB]->SetLineStyle(kDashed);

  canv->Divide(2,1);

  canv->cd(1);
  Int_t f = dist[qA]->GetMaximum() > dist[qB]->GetMaximum() ? qA:qB;
  dist[f]->Draw();
  dist[(f+1)%2]->Draw("SAME");

  canv->cd(2);
  corr->Draw("colz");
  canv->GetPad(2)->SetGrid(1,1);
};

