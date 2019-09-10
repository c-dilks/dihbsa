R__LOAD_LIBRARY(DihBsa)

#include "Constants.h"
#include "EventTree.h"
#include "Tools.h"
#include "TString.h"
#include "TMath.h"

enum fenum {kGen,kRec};
int f;
EventTree * ev[2];

void MCmatch(TString fname=
  "out_clasdispr.00.e10.600.emn0.75tmn.09.xs80.53nb.dis.0062.nrad.dat.evio.hipo.root") {

  Int_t pairType = EncodePairType(kPip,kPim); // make sure ordering obeys dihHadIdx...

  TString infileN[2];
  infileN[kGen] = "outroot.MC.gen/";
  infileN[kRec] = "outroot.MC.rec/";
  for(f=0; f<2; f++) {
    infileN[f] += fname;
    ev[f] = new EventTree(infileN[f],pairType);
  };

  Bool_t success = ev[kRec]->BuildEvnumMap();
  if(!success) return;

  for(int i=0; i<ev[kGen]->ENT; i++) {
    ev[kGen]->GetEvent(i);
    if(pairType==ev[kGen]->pairType) {
      if(ev[kRec]->FindEvent(ev[kGen]->evnum)) {
        printf("event %d\n",ev[kGen]->evnum);
        printf("GEN: hadE[qA]=%.3f  hadE[qB]=%.3f\n",
          ev[kGen]->hadE[qA],ev[kGen]->hadE[qB]);
        printf("REC: hadE[qA]=%.3f  hadE[qB]=%.3f\n",
          ev[kRec]->hadE[qA],ev[kRec]->hadE[qB]);
      };
    };
  };
};
