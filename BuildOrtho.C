R__LOAD_LIBRARY(DihBsa)

#include "Constants.h"
#include "EventTree.h"

void BuildOrtho(TString inDir="outroot.fall18.someMore") {
  const Int_t NBINS = 50;
  TFile * outfile = new TFile("ortho.root","RECREATE");

  // yield distributions
  TH2D * d2 = new TH2D("d2","d2",NBINS,-PI,PI,NBINS,-PI,PI);
  TH3D * d3 = new TH3D("d3","d3",NBINS,-PI,PI,NBINS,-PI,PI,NBINS,0,PI);

  // test plots
  TH2D * pph[2];
  pph[0] = new TH2D("pphA","#phi_{h}^{#pi1#pi2} vs. #phi_{h}^{#pi1}",NBINS,-PI,PI,
                                                                     NBINS,-PI,PI);
  pph[1] = new TH2D("pphB","#phi_{h}^{#pi1#pi2} vs. #phi_{h}^{#pi2}",NBINS,-PI,PI,
                                                                     NBINS,-PI,PI);

  EventTree * ev = new EventTree(TString(inDir+"/*.root"),EncodePairType(kPip,kPim));
  for(int i=0; i<ev->ENT; i++) {
    ev->GetEvent(i);
    if(ev->Valid()) {
      d2->Fill(ev->PhiR,ev->PhiH);
      d3->Fill(ev->PhiR,ev->PhiH,ev->theta);
      for(int h=0; h<2; h++) {
        pph[h]->Fill(ev->GetDihadronObj()->GetSingleHadronPhiH(h),ev->PhiH);
      };
    };
  };
  d2->Write();
  d3->Write();
  for(int h=0; h<2; h++) pph[h]->Write();
};
