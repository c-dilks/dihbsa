void BuildOrthoStefan(TString infileN="epipX.root") {

  TFile * infile = new TFile(infileN,"READ");
  TTree * tr = (TTree*) infile->Get("events");

  const Int_t NBINS = 50;
  TFile * outfile = new TFile("ortho.root","RECREATE");
  TH1D * d1 = new TH1D("d1","d1",
    NBINS,-TMath::Pi(),TMath::Pi());

  tr->Project("d1","phi","");
  //tr->Project("d1","phi","pt>1");
  //tr->Project("d1","phi","pt<1");

  d1->Write();
};
