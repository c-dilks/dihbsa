// plot ratio of specific plots in eicPlots*.root files

void RatioPlot(
  TString numerFileN="eicPlots.pythia_5x41_smear.ymin_0.03.root",
  TString denomFileN="eicPlots.pythia_5x41_smear.ymin_0.00.root",
  TString particle="hadronB",
  TString outfileN="PperpRatio_PiMinus_ymin0.03"
) {

  enum numerdenom {num,den};
  TFile * infile[2];
  infile[num] = new TFile(numerFileN,"READ");
  infile[den] = new TFile(denomFileN,"READ");

  const Int_t nBins = 9;
  TH1D * dist[2][nBins];
  TH1D * ratio[nBins];
  TString distN;
  TCanvas * canv[nBins];
  for(int b=0; b<nBins; b++) {
    for(int f=0; f<2; f++) {
      distN = "singlePlots/"+particle+
        "_PperpDistLin_"+Form("%d",b);
      dist[f][b] = (TH1D*) infile[f]->Get(distN);
      printf("f%d b%d %s @ %p\n",f,b,
        distN.Data(),(void*)dist[f][b]);
      dist[f][b]->Sumw2();
    };
    ratio[b] = (TH1D*) dist[num][b]->Clone();
    ratio[b]->Divide(dist[den][b]);
    canv[b] = new TCanvas(Form("canv%d",b),Form("canv%d",b),
      1600,700);
    canv[b]->Divide(2,1);
    canv[b]->cd(1);
    canv[b]->GetPad(1)->SetGrid(1,1);
    dist[num][b]->SetLineColor(kBlue);
    dist[den][b]->SetLineColor(kRed);
    dist[den][b]->Draw("ERR");
    dist[num][b]->Draw("ERR SAME");
    canv[b]->cd(2);
    canv[b]->GetPad(2)->SetGrid(1,1);
    ratio[b]->SetLineColor(kBlack);
    ratio[b]->Draw("ERR");
    canv[b]->Print(Form("%s_bin%d.pdf",outfileN.Data(),b),"pdf");
  };

  //for(int f=0; f<2; f++) infile[f]->Close();
};
