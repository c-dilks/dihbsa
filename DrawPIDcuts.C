void DrawPIDcuts(TString infileN = "eicPlots.pythia_18x275_smear.root") {
  TFile * infile = new TFile(infileN,"READ");
  TH2D * dist = (TH2D*) infile->Get("singlePlots/hadronA_EtaVsP_0");
  TCanvas * canv = new TCanvas("canv","canv",800,700);
  canv->SetLogx();
  //canv->SetLogz();
  canv->SetGrid(1,1);
  dist->Draw("colz");
  
  const int N=5;
  TLine * line[N];
  line[0] = new TLine(7,-3.5,7,-1);
  line[1] = new TLine(5,-1,5,1);
  line[2] = new TLine(8,1,8,2);
  line[3] = new TLine(20,2,20,3);
  line[4] = new TLine(45,3,45,3.5);
  for(int i=0; i<N; i++) {
    line[i]->SetLineColor(kBlack);
    line[i]->SetLineWidth(10);
    line[i]->Draw();
  };
};
