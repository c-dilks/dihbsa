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
#include "EventTree.h"

Float_t ModAngle(Float_t ang);
void HadronCompareCanv(TCanvas * canv, TH1F * dist[2], TH2F * corr);

int main(int argc, char** argv) {

   // ARGUMENTS
   TString inDir = "outroot";
   if(argc>1) inDir = TString(argv[1]);

   gSystem->Load("src/DihBsa.so");
   EventTree * ev = new EventTree(TString(inDir+"/out*.root"));


   TFile * outfile = new TFile("phiRplots.root","RECREATE");


   const Int_t NBINS = 100;
   Float_t distMax = 1.1;

   TH1F * dist_T_byKt = new TH1F("dist_T_byKt",
     "sin[ #phi_{R}(T,k_{T}) ] distribution -- using R_{T} via k_{T} formula;sin[ #phi_{R} ]",
     NBINS,-distMax,distMax);
   TH1F * dist_T_byRej = new TH1F("dist_T_byRej",
     "sin[ #phi_{R}(T,rej) ] distribution -- using R_{T} via vector rejection;sin[ #phi_{R} ]",
     NBINS,-distMax,distMax);
   TH1F * dist_Perp = new TH1F("dist_Perp",
     "sin[ #phi_{R}(#perp  ,rej) ] distribution -- using R_{ #perp  } via vector rejection;sin[ #phi_{R} ]",
     NBINS,-distMax,distMax);
   dist_T_byKt->SetLineColor(kRed);
   dist_T_byRej->SetLineColor(kGreen+1);
   dist_Perp->SetLineColor(kBlue);

   TH1F * diff__T_byKt__T_byRej = new TH1F("diff__T_byKt__T_byRej",
     "sin[ #phi_{R}(T,k_{T}) - #phi_{R}(T,rej) ]",
     NBINS,-distMax,distMax);
   TH1F * diff__T_byKt__Perp = new TH1F("diff__T_byKt__Perp",
     "sin[ #phi_{R}(T,k_{T}) - #phi_{R}(#perp  ,rej) ]",
     NBINS,-distMax,distMax);
   TH1F * diff__T_byRej__Perp = new TH1F("diff__T_byRej__Perp",
     "sin[ #phi_{R}(T,rej) - #phi_{R}(#perp  ,rej) ]",
     NBINS,-distMax,distMax);
   diff__T_byKt__T_byRej->SetLineColor(kRed);
   diff__T_byKt__Perp->SetLineColor(kGreen+1);
   diff__T_byRej__Perp->SetLineColor(kBlue);

   TH2F * corr__T_byKt__T_byRej = new TH2F("corr__T_byKt__T_byRej",
     "sin[ #phi_{R}(T,k_{T}) ] vs. sin[ #phi_{R}(T,rej) ]",
     NBINS,-distMax,distMax,
     NBINS,-distMax,distMax);
   TH2F * corr__T_byKt__Perp = new TH2F("corr__T_byKt__Perp",
     "sin[ #phi_{R}(T,k_{T}) ] vs. sin[ #phi_{R}(#perp  ,rej) ]",
     NBINS,-distMax,distMax,
     NBINS,-distMax,distMax);
   TH2F * corr__T_byRej__Perp = new TH2F("corr__T_byRej__Perp",
     "sin[ #phi_{R}(T,rej) ] vs. sin[ #phi_{R}(#perp  ,rej) ]",
     NBINS,-distMax,distMax,
     NBINS,-distMax,distMax);

   TH2F * breitVsLab_T_byKt = new TH2F("breitVsLab_T_byKt",
     "sin[ #phi_{R}(T,k_{T}) ] via Breit Frame vs. via Lab Frame",
     NBINS,-distMax,distMax,
     NBINS,-distMax,distMax);
   TH2F * breitVsLab_T_byRej = new TH2F("breitVsLab_T_byRej",
     "sin[ #phi_{R}(T,Rej) ] via Breit Frame vs. via Lab Frame",
     NBINS,-distMax,distMax,
     NBINS,-distMax,distMax);
   TH2F * breitVsLab_Perp = new TH2F("breitVsLab_Perp",
     "sin[ #phi_{R}(#perp  ,rej) ] via Breit Frame vs. via Lab Frame",
     NBINS,-distMax,distMax,
     NBINS,-distMax,distMax);


   Float_t ang_T_byKt;
   Float_t ang_T_byRej;
   Float_t ang_Perp;
   Float_t angDiff;

   Float_t b_ang_T_byKt;
   Float_t b_ang_T_byRej;
   Float_t b_ang_Perp;

   printf("begin loop through %d events...\n",ev->ENT);
   for(int i=0; i<ev->ENT; i++) {

     ev->GetEvent(i);

     if(ev->cutDihadron && ev->cutQ2 && ev->cutW && ev->cutY) {

       // lab frame
       ang_T_byKt = TMath::Sin(ev->PhiR_T_byKt);
       ang_T_byRej = TMath::Sin(ev->PhiR_T_byRej);
       ang_Perp = TMath::Sin(ev->PhiR_Perp);

       // breit frame
       b_ang_T_byKt = TMath::Sin(ev->b_PhiR_T_byKt);
       b_ang_T_byRej = TMath::Sin(ev->b_PhiR_T_byRej);
       b_ang_Perp = TMath::Sin(ev->b_PhiR_Perp);


       dist_T_byKt->Fill(ang_T_byKt);
       dist_T_byRej->Fill(ang_T_byRej);
       dist_Perp->Fill(ang_Perp);


       corr__T_byKt__T_byRej->Fill(ang_T_byRej,ang_T_byKt);
       corr__T_byKt__Perp->Fill(ang_Perp,ang_T_byKt);
       corr__T_byRej__Perp->Fill(ang_Perp,ang_T_byRej);


       angDiff = TMath::Sin(ModAngle(ang_T_byKt - ang_T_byRej));
       diff__T_byKt__T_byRej->Fill(angDiff);

       angDiff = TMath::Sin(ModAngle(ang_T_byKt - ang_Perp));
       diff__T_byKt__Perp->Fill(angDiff);

       angDiff = TMath::Sin(ModAngle(ang_T_byRej - ang_Perp));
       diff__T_byRej__Perp->Fill(angDiff);


       breitVsLab_T_byKt->Fill(ang_T_byKt,b_ang_T_byKt);
       breitVsLab_T_byRej->Fill(ang_T_byRej,b_ang_T_byRej);
       breitVsLab_Perp->Fill(ang_Perp,b_ang_Perp);
     };


   };

   dist_T_byKt->Write();
   dist_T_byRej->Write();
   dist_Perp->Write();

   corr__T_byKt__T_byRej->Write();
   corr__T_byKt__Perp->Write();
   corr__T_byRej__Perp->Write();

   diff__T_byKt__T_byRej->Write();
   diff__T_byKt__Perp->Write();
   diff__T_byRej__Perp->Write();

   breitVsLab_T_byKt->Write();
   breitVsLab_T_byRej->Write();
   breitVsLab_Perp->Write();


   outfile->Close();

};


Float_t ModAngle(Float_t ang) {
  while(ang>PI) ang-=2*PI;
  while(ang<-PI) ang+=2*PI;
  return ang;
};
