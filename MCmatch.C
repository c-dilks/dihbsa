R__LOAD_LIBRARY(DihBsa)

#include "Constants.h"
#include "EventTree.h"
#include "Tools.h"
#include "TString.h"
#include "TMath.h"
#include "TTree.h"

enum fenum {kGen,kRec};
int f;
EventTree * ev[2];

void MCmatch(TString fname=
  "out_clasdispr.00.e10.600.emn0.75tmn.09.xs80.53nb.dis.0410.nrad.dat.evio.hipo.root"
) {

  Int_t pairType = EncodePairType(kPip,kPim); // make sure ordering obeys dihHadIdx...
  
  TFile * outfile = new TFile(TString("match."+fname),"RECREATE");

  TString mgStr[2];
  mgStr[kGen] = "gen";
  mgStr[kRec] = "rec";


  TString infileN[2];
  for(f=0; f<2; f++) {
    infileN[f] = "outroot.MC."+mgStr[f]+"/"+fname;
    ev[f] = new EventTree(infileN[f],pairType);
  };

  Float_t diff_hadE[2];
  Float_t diff_hadP[2];
  Float_t diff_PhiH;
  Float_t diff_PhiR;

  TTree * mtr = new TTree("mtr","mtr");
  mtr->Branch("evnum",&(ev[kRec]->evnum),"evnum/I");
  for(f=0; f<2; f++) {
    mtr->Branch(TString(mgStr[f]+"_hadE"), ev[f]->hadE,
                TString(mgStr[f]+"_hadE[2]/F"));
    mtr->Branch(TString(mgStr[f]+"_hadP"), ev[f]->hadP,
                TString(mgStr[f]+"_hadP[2]/F"));
    mtr->Branch(TString(mgStr[f]+"_helicity"), &(ev[f]->helicity),
                TString(mgStr[f]+"_helicity/I"));
    mtr->Branch(TString(mgStr[f]+"_hadOrder"), &(ev[f]->hadOrder),
                TString(mgStr[f]+"_hadOrder/I"));
  };
  mtr->Branch("diff_hadP",diff_hadP,"diff_hadP[2]/F");
  mtr->Branch("diff_PhiH",&diff_PhiH,"diff_PhiH/F");
  mtr->Branch("diff_PhiR",&diff_PhiR,"diff_PhiR/F");

  Bool_t success = ev[kGen]->BuildMatchTable();
  if(!success) return;

  Int_t nTotal=0;
  Int_t nMatches=0;


  for(int i=0; i<ev[kRec]->ENT; i++) {
    ev[kRec]->GetEvent(i);
    if(pairType==ev[kRec]->pairType) {

      nTotal++;

      if( ev[kGen]->FindEvent( ev[kRec]->evnum,
                               ev[kRec]->hadP[qA],
                               ev[kRec]->hadP[qB]) ) {
        nMatches++;

        for(int h=0; h<2; h++) {
          diff_hadE[h] = (ev[kGen]->hadE[h] - ev[kRec]->hadE[h]) / ev[kRec]->hadE[h]; 
          diff_hadP[h] = (ev[kGen]->hadP[h] - ev[kRec]->hadP[h]) / ev[kRec]->hadP[h]; 
        };
        diff_PhiH = Tools::AdjAngle(ev[kGen]->PhiH - ev[kRec]->PhiH) / ev[kRec]->PhiH;
        diff_PhiR = Tools::AdjAngle(ev[kGen]->PhiR - ev[kRec]->PhiR) / ev[kRec]->PhiR;

        mtr->Fill();

      };
    };
  };
  mtr->Write();

  new TCanvas();
  Int_t nbins = 40;
  Float_t drawLim = 0.1;
  TString drawStr;
  drawStr = Form(">>d(%d,%f,%f,%d,%f,%f)",nbins,-drawLim,drawLim,nbins,-drawLim,drawLim);
  mtr->Draw(TString("diff_hadP[0]:diff_hadP[1]"+drawStr),"","colz");

  new TCanvas();
  drawLim = 1;
  drawStr = Form(">>e(%d,%f,%f,%d,%f,%f)",nbins,-drawLim,drawLim,nbins,-drawLim,drawLim);
  mtr->Draw(TString("diff_PhiH:diff_PhiR"+drawStr),"","colz");

  printf("\n\n%d / %d dihadron pairs match (%.3f%%)\n\n\n",nMatches,nTotal,
    Float_t(100*nMatches)/nTotal);
};
