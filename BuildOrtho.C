R__LOAD_LIBRARY(DihBsa)

#include "Constants.h"
#include "EventTree.h"

void BuildOrtho(TString inDir="outroot.fall18.some") {
  const Int_t NBINS = 50;
  TFile * outfile = new TFile("ortho.root","RECREATE");
  TH2D * d2 = new TH2D("d2","d2",NBINS,-PI,PI,NBINS,-PI,PI);
  TH3D * d3 = new TH3D("d3","d3",NBINS,-PI,PI,NBINS,-PI,PI,NBINS,0,PI);
  EventTree * ev = new EventTree(TString(inDir+"/*.root"),EncodePairType(kPip,kPim));
  for(int i=0; i<ev->ENT; i++) {
    ev->GetEvent(i);
    if(ev->Valid()) {
      d2->Fill(ev->PhiR,ev->PhiH);
      d3->Fill(ev->PhiR,ev->PhiH,ev->theta);
    };
  };
  d2->Write();
  d3->Write();
};
