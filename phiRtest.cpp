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
   Bool_t useBreit = false;
   if(argc>1) inDir = TString(argv[1]);
   if(argc>2) useBreit = (Bool_t)strtof(argv[2],NULL);

   gSystem->Load("src/DihBsa.so");
   EventTree * ev = new EventTree(TString(inDir+"/out*.root"));
   TString frame = useBreit ? " --- Breit frame":" --- Lab frame ";
   printf("\nPhiR definition comparison plots will be in %s\n",frame.Data());


   TFile * outfile = new TFile("phiRplots.root","RECREATE");


   const Int_t NBINS = 100;
   Float_t distMax = 1.1;

   TH1F * distPhiRp = new TH1F("distPhiRp",
     "sin[ #phi_{R}(T,k_{T}) ] distribution -- using R_{T} via k_{T} formula;sin[ #phi_{R} ]",
     NBINS,-distMax,distMax);
   TH1F * distPhiRp_r = new TH1F("distPhiRp_r",
     "sin[ #phi_{R}(T,rej) ] distribution -- using R_{T} via vector rejection;sin[ #phi_{R} ]",
     NBINS,-distMax,distMax);
   TH1F * distPhiRq = new TH1F("distPhiRq",
     "sin[ #phi_{R}(#perp  ,rej) ] distribution -- using R_{ #perp  } via vector rejection;sin[ #phi_{R} ]",
     NBINS,-distMax,distMax);
   distPhiRp->SetLineColor(kRed);
   distPhiRp_r->SetLineColor(kGreen+1);
   distPhiRq->SetLineColor(kBlue);

   TH1F * diff_PhiRp_PhiRp_r = new TH1F("diff_PhiRp_PhiRp_r",
     "sin[ #phi_{R}(T,k_{T}) - #phi_{R}(T,rej) ]",
     NBINS,-distMax,distMax);
   TH1F * diff_PhiRp_PhiRq = new TH1F("diff_PhiRp_PhiRq",
     "sin[ #phi_{R}(T,k_{T}) - #phi_{R}(#perp  ,rej) ]",
     NBINS,-distMax,distMax);
   TH1F * diff_PhiRp_r_PhiRq = new TH1F("diff_PhiRp_r_PhiRq",
     "sin[ #phi_{R}(T,rej) - #phi_{R}(#perp  ,rej) ]",
     NBINS,-distMax,distMax);
   diff_PhiRp_PhiRp_r->SetLineColor(kRed);
   diff_PhiRp_PhiRq->SetLineColor(kGreen+1);
   diff_PhiRp_r_PhiRq->SetLineColor(kBlue);

   TH2F * corr_PhiRp_PhiRp_r = new TH2F("corr_PhiRp_PhiRp_r",
     "sin[ #phi_{R}(T,k_{T}) ] vs. sin[ #phi_{R}(T,rej) ]",
     NBINS,-distMax,distMax,
     NBINS,-distMax,distMax);
   TH2F * corr_PhiRp_PhiRq = new TH2F("corr_PhiRp_PhiRq",
     "sin[ #phi_{R}(T,k_{T}) ] vs. sin[ #phi_{R}(#perp  ,rej) ]",
     NBINS,-distMax,distMax,
     NBINS,-distMax,distMax);
   TH2F * corr_PhiRp_r_PhiRq = new TH2F("corr_PhiRp_r_PhiRq",
     "sin[ #phi_{R}(T,rej) ] vs. sin[ #phi_{R}(#perp  ,rej) ]",
     NBINS,-distMax,distMax,
     NBINS,-distMax,distMax);




   TString newTitle;

   newTitle = distPhiRp->GetTitle() + frame;
   distPhiRp->SetTitle(newTitle);

   newTitle = distPhiRp_r->GetTitle() + frame;
   distPhiRp_r->SetTitle(newTitle);

   newTitle = distPhiRq->GetTitle() + frame;
   distPhiRq->SetTitle(newTitle);

   newTitle = diff_PhiRp_PhiRp_r->GetTitle() + frame;
   diff_PhiRp_PhiRp_r->SetTitle(newTitle);

   newTitle = diff_PhiRp_PhiRq->GetTitle() + frame;
   diff_PhiRp_PhiRq->SetTitle(newTitle);

   newTitle = diff_PhiRp_r_PhiRq->GetTitle() + frame;
   diff_PhiRp_r_PhiRq->SetTitle(newTitle);

   newTitle = corr_PhiRp_PhiRp_r->GetTitle() + frame;
   corr_PhiRp_PhiRp_r->SetTitle(newTitle);

   newTitle = corr_PhiRp_PhiRq->GetTitle() + frame;
   corr_PhiRp_PhiRq->SetTitle(newTitle);

   newTitle = corr_PhiRp_r_PhiRq->GetTitle() + frame;
   corr_PhiRp_r_PhiRq->SetTitle(newTitle);



   TH2F * breitVsLabPhiRp = new TH2F("breitVsLabPhiRp",
     "sin[ #phi_{R}(T,k_{T}) ] via Breit Frame vs. via Lab Frame",
     NBINS,-distMax,distMax,
     NBINS,-distMax,distMax);
   TH2F * breitVsLabPhiRp_r = new TH2F("breitVsLabPhiRp_r",
     "sin[ #phi_{R}(T,Rej) ] via Breit Frame vs. via Lab Frame",
     NBINS,-distMax,distMax,
     NBINS,-distMax,distMax);
   TH2F * breitVsLabPhiRq = new TH2F("breitVsLabPhiRq",
     "sin[ #phi_{R}(#perp  ,rej) ] via Breit Frame vs. via Lab Frame",
     NBINS,-distMax,distMax,
     NBINS,-distMax,distMax);


   TString plotTitle;
   TString sigmaStr = "sin[ #phi_{R}(T,k_{T}) - #phi_{R}(#perp  ,rej) ]";
   plotTitle = sigmaStr+" vs. Q^{~2}";
   TH2F * sigmaVsQ2 = new TH2F("sigmaVsQ2",plotTitle,
     NBINS,0,10,
     NBINS,-distMax,distMax);
   plotTitle = sigmaStr+" vs. P_{ h}";
   TH2F * sigmaVsPh = new TH2F("sigmaVsPh",plotTitle,
     NBINS,0,8,
     NBINS,-distMax,distMax);
   plotTitle = sigmaStr+" vs. P_{ h#perp}";
   TH2F * sigmaVsPhPerp = new TH2F("sigmaVsPhPerp",plotTitle,
     NBINS,0,2,
     NBINS,-distMax,distMax);
   plotTitle = sigmaStr+" vs. z";
   TH2F * sigmaVsZ = new TH2F("sigmaVsZ",plotTitle,
     NBINS,0,1,
     NBINS,-distMax,distMax);
   plotTitle = sigmaStr+" vs. x_{F}";
   TH2F * sigmaVsXF = new TH2F("sigmaVsXF",plotTitle,
     NBINS*2,-1,1,
     NBINS,-distMax,distMax);
   TH1F * sigmaDist = new TH1F("sigmaDist",sigmaStr,
     NBINS,-distMax,distMax);





   Float_t angPhiRp;
   Float_t angPhiRp_r;
   Float_t angPhiRq;

   Float_t l_angPhiRp;
   Float_t l_angPhiRp_r;
   Float_t l_angPhiRq;

   Float_t b_angPhiRp;
   Float_t b_angPhiRp_r;
   Float_t b_angPhiRq;

   Float_t angDiff;
   Float_t sigma;

   printf("begin loop through %d events...\n",ev->ENT);
   for(int i=0; i<ev->ENT; i++) {

     ev->GetEvent(i);

     if(ev->cutDihadron && ev->cutQ2 && ev->cutW && ev->cutY) {

       // lab frame
       l_angPhiRp = ev->PhiRp;
       l_angPhiRp_r = ev->PhiRp_r;
       l_angPhiRq = ev->PhiRq;

       // breit frame
       b_angPhiRp = ev->b_PhiRp;
       b_angPhiRp_r = ev->b_PhiRp_r;
       b_angPhiRq = ev->b_PhiRq;

       if(useBreit) {
         angPhiRp = b_angPhiRp;
         angPhiRp_r = b_angPhiRp_r;
         angPhiRq = b_angPhiRq;
       } else {
         angPhiRp = l_angPhiRp;
         angPhiRp_r = l_angPhiRp_r;
         angPhiRq = l_angPhiRq;
       };


       distPhiRp->Fill(TMath::Sin(angPhiRp));
       distPhiRp_r->Fill(TMath::Sin(angPhiRp_r));
       distPhiRq->Fill(TMath::Sin(angPhiRq));


       corr_PhiRp_PhiRp_r->Fill(TMath::Sin(angPhiRp_r),TMath::Sin(angPhiRp));
       corr_PhiRp_PhiRq->Fill(TMath::Sin(angPhiRq),TMath::Sin(angPhiRp));
       corr_PhiRp_r_PhiRq->Fill(TMath::Sin(angPhiRq),TMath::Sin(angPhiRp_r));


       angDiff = TMath::Sin(ModAngle(angPhiRp - angPhiRp_r));
       diff_PhiRp_PhiRp_r->Fill(angDiff);

       angDiff = TMath::Sin(ModAngle(angPhiRp - angPhiRq));
       sigma = angDiff;
       diff_PhiRp_PhiRq->Fill(angDiff);

       angDiff = TMath::Sin(ModAngle(angPhiRp_r - angPhiRq));
       diff_PhiRp_r_PhiRq->Fill(angDiff);


       breitVsLabPhiRp->Fill(TMath::Sin(l_angPhiRp),TMath::Sin(b_angPhiRp));
       breitVsLabPhiRp_r->Fill(TMath::Sin(l_angPhiRp_r),TMath::Sin(b_angPhiRp_r));
       breitVsLabPhiRq->Fill(TMath::Sin(l_angPhiRq),TMath::Sin(b_angPhiRq));

       sigmaVsQ2->Fill(ev->Q2,sigma);
       sigmaVsPh->Fill(ev->Ph,sigma);
       sigmaVsPhPerp->Fill(ev->PhPerp,sigma);
       sigmaVsZ->Fill(ev->Zpair,sigma);
       sigmaVsXF->Fill(ev->xF,sigma);
       if(ev->Q2 > 2 && ev->Q2 < 3) sigmaDist->Fill(sigma);
     };


   };

   distPhiRp->Write();
   distPhiRp_r->Write();
   distPhiRq->Write();

   corr_PhiRp_PhiRp_r->Write();
   corr_PhiRp_PhiRq->Write();
   corr_PhiRp_r_PhiRq->Write();

   diff_PhiRp_PhiRp_r->Write();
   diff_PhiRp_PhiRq->Write();
   diff_PhiRp_r_PhiRq->Write();

   breitVsLabPhiRp->Write();
   breitVsLabPhiRp_r->Write();
   breitVsLabPhiRq->Write();

   sigmaVsQ2->Write();
   sigmaVsPh->Write();
   sigmaVsPhPerp->Write();
   sigmaVsZ->Write();
   sigmaVsXF->Write();
   sigmaDist->Write();


   outfile->Close();

};


Float_t ModAngle(Float_t ang) {
  while(ang>PI) ang-=2*PI;
  while(ang<-PI) ang+=2*PI;
  return ang;
};
