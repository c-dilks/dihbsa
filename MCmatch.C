// plots asymmetry modulation versus various other kinematics, to get an overall sense
// of the acceptance from the data

R__LOAD_LIBRARY(DihBsa)

#include "Constants.h"
#include "EventTree.h"
#include "Dihadron.h"
#include "Tools.h"
#include "TString.h"
#include "TMath.h"
#include "TTree.h"
#include "TH1.h"

enum fenum {kGen,kRec};
int f;
EventTree * ev[2];

void MCmatch(TString fname=
"clasdispr.00.e10.600.emn0.75tmn.09.xs81.61nb.dis.0410.dat.hipo.root"
) {

  Int_t pairType = EncodePairType(kPip,kPim); // make sure ordering obeys dihHadIdx...
  
  TFile * outfile = new TFile(TString("match."+fname),"RECREATE");

  TString mgStr[2];
  mgStr[kGen] = "gen";
  mgStr[kRec] = "rec";

  // instantiate EventTree
  TString infileN[2];
  for(f=0; f<2; f++) {
    infileN[f] = "outroot.MC."+mgStr[f]+"/"+fname;
    ev[f] = new EventTree(infileN[f],pairType);
  };

  Float_t delta_hadE[2];
  Float_t delta_hadP[2];
  Float_t diff_PhiH;
  Float_t diff_PhiR;

  // define output tree 'mtr'
  TTree * mtr = new TTree("mtr","mtr");
  mtr->Branch("evnum",&(ev[kRec]->evnum),"evnum/I");
  for(f=0; f<2; f++) {
    mtr->Branch(TString(mgStr[f]+"_hadE"), ev[f]->hadE,
                TString(mgStr[f]+"_hadE[2]/F"));
    mtr->Branch(TString(mgStr[f]+"_hadP"), ev[f]->hadP,
                TString(mgStr[f]+"_hadP[2]/F"));
    mtr->Branch(TString(mgStr[f]+"_hadOrder"), &(ev[f]->hadOrder),
                TString(mgStr[f]+"_hadOrder/I"));
  };
  mtr->Branch("delta_hadP",delta_hadP,"delta_hadP[2]/F");
  mtr->Branch("delta_hadE",delta_hadE,"delta_hadE[2]/F");
  mtr->Branch("diff_PhiH",&diff_PhiH,"diff_PhiH/F");
  mtr->Branch("diff_PhiR",&diff_PhiR,"diff_PhiR/F");
  mtr->Branch("matchDiff",&(ev[kGen]->matchDiff),"matchDiff/F");

  // define matching fraction histos
  // -- in a single execution of MCmatch.C, these histograms are just distributions
  //    of the kinematics; they aren't yet the matching fractions
  // -- to get matching fractions, call loopMCmatch.sh, which then calls
  //    drawMatchFraction.C, which draws the plots
  TH1D * MhMF[2];
  TH1D * XMF[2];
  TH1D * ZMF[2];
  TString mStr[2];
  mStr[0] = "All";
  mStr[1] = "Matched";
  Int_t NBINS = 7;
  for(int m=0; m<2; m++) {
    MhMF[m] = new TH1D(TString("MhMF"+mStr[m]),"Matching Fraction vs. M_{h}",NBINS,0,3);
    XMF[m] = new TH1D(TString("XMF"+mStr[m]),"Matching Fraction vs. x",NBINS,0,1);
    ZMF[m] = new TH1D(TString("ZMF"+mStr[m]),"Matching Fraction vs. z_{pair}",NBINS,0,1);
    MhMF[m]->Sumw2();
    XMF[m]->Sumw2();
    ZMF[m]->Sumw2();
  };


  // build match table of MCgen events
  Bool_t success = ev[kGen]->BuildMatchTable();
  if(!success) return;
  Bool_t matchFound;

  Int_t nTotal=0;
  Int_t nMatches=0;

  Float_t matchDiff;
  TH1D * matchDiffDist = new TH1D("matchDiffDist",
    "D distribution (no D cut)",300,0,2);
  TH2D * matchDiffVsPdiff = new TH2D("matchDiffVsPdiff",
    "D vs. #Delta P_{h} (no D cut);#Delta P_{h};D",300,-0.3,0.3,300,0,2);


  // loop through MCrec events
  for(int i=0; i<ev[kRec]->ENT; i++) {
    ev[kRec]->GetEvent(i);

    // check dihadron cuts
    if(pairType==ev[kRec]->pairType) {
    //if(pairType==ev[kRec]->pairType && ev[kRec]->Valid()) {

      nTotal++;

      // look for matching MCgen dihadron
      matchFound = ev[kGen]->FindEvent( ev[kRec]->evnum, ev[kRec]->GetDihadronObj() );
      if(ev[kGen]->matchDiff<1e6) {
        matchDiffDist->Fill(ev[kGen]->matchDiff);
        matchDiffVsPdiff->Fill( 
          (ev[kGen]->Ph - ev[kRec]->Ph) / ev[kRec]->Ph, 
          ev[kGen]->matchDiff );
      };

      if(matchFound) {
        nMatches++;

        // fill plots and mtr
        for(int h=0; h<2; h++) {
          delta_hadE[h] = (ev[kGen]->hadE[h] - ev[kRec]->hadE[h]) / ev[kRec]->hadE[h]; 
          delta_hadP[h] = (ev[kGen]->hadP[h] - ev[kRec]->hadP[h]) / ev[kRec]->hadP[h]; 
        };
        diff_PhiH = Tools::AdjAngle(ev[kGen]->PhiH - ev[kRec]->PhiH);
        diff_PhiR = Tools::AdjAngle(ev[kGen]->PhiR - ev[kRec]->PhiR);

        mtr->Fill();

        MhMF[1]->Fill(ev[kRec]->Mh); // fill "matched" distribution (MF numerator)
        XMF[1]->Fill(ev[kRec]->x);
        ZMF[1]->Fill(ev[kRec]->Zpair);

      };

      MhMF[0]->Fill(ev[kRec]->Mh); // fill "all" distribution (MF denominator)
      XMF[0]->Fill(ev[kRec]->x);
      ZMF[0]->Fill(ev[kRec]->Zpair);

    };
  };


  //MhMF[1]->Divide(MhMF[0]);
  //XMF[1]->Divide(XMF[0]);
  //ZMF[1]->Divide(ZMF[0]);
  for(int m=0; m<2; m++) MhMF[m]->Write();
  for(int m=0; m<2; m++) XMF[m]->Write();
  for(int m=0; m<2; m++) ZMF[m]->Write();
  matchDiffDist->Write();
  matchDiffVsPdiff->Write();

  mtr->Write();

  new TCanvas();
  Int_t nbins = 40;
  Float_t drawLim = 0.1;
  TString drawStr;
  drawStr = Form(">>d(%d,%f,%f,%d,%f,%f)",nbins,-drawLim,drawLim,nbins,-drawLim,drawLim);
  mtr->Draw(TString("delta_hadP[0]:delta_hadP[1]"+drawStr),"","colz");

  new TCanvas();
  drawLim = 1;
  drawStr = Form(">>e(%d,%f,%f,%d,%f,%f)",nbins,-drawLim,drawLim,nbins,-drawLim,drawLim);
  mtr->Draw(TString("diff_PhiH:diff_PhiR"+drawStr),"","colz");

  printf("\n\n%d / %d dihadron pairs match (%.3f%%)\n\n\n",nMatches,nTotal,
    Float_t(100*nMatches)/nTotal);
};
