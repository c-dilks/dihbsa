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

  TTree * mtr = new TTree("mtr","mtr");
  mtr->Branch("evnum",&(ev[kGen]->evnum),"evnum/I");
  for(f=0; f<2; f++) {
    mtr->Branch(TString(mgStr[f]+"_hadE"), ev[f]->hadE,
                TString(mgStr[f]+"_hadE[2]/F"));
    mtr->Branch(TString(mgStr[f]+"_helicity"), &(ev[f]->helicity),
                TString(mgStr[f]+"_helicity/I"));
  };
  mtr->Branch("diff_hadE",diff_hadE,"diff_hadE[2]/F");

  Bool_t success = ev[kRec]->BuildEvnumMap();
  if(!success) return;


  for(int i=0; i<ev[kGen]->ENT; i++) {
    ev[kGen]->GetEvent(i);
    if(pairType==ev[kGen]->pairType) {
      if(ev[kRec]->FindEvent(ev[kGen]->evnum)) {
        for(int h=0; h<2; h++) {
          diff_hadE[h] = 
            TMath::Abs(ev[kGen]->hadE[h] - ev[kRec]->hadE[h]) / ev[kGen]->hadE[h]; 
        };

        /*
        printf("event %d\n",ev[kGen]->evnum);
        printf("GEN: hadE[qA]=%.3f  hadE[qB]=%.3f\n",
          ev[kGen]->hadE[qA],ev[kGen]->hadE[qB]);
        printf("REC: hadE[qA]=%.3f  hadE[qB]=%.3f\n",
          ev[kRec]->hadE[qA],ev[kRec]->hadE[qB]);
        printf("DIFF: diff_hadE[qA]=%.3f  diff_hadE[qB]=%.3f\n",
          diff_hadE[qA],diff_hadE[qB]);
        */

        mtr->Fill();

      };
    };
  };
  mtr->Write();
};
