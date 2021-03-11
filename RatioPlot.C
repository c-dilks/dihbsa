// plot ratio of specific plots in eicPlots*.root files
//
// WARNING: drawing histograms with variable bin widths (e.g.,
// even bins in log scale), does not work if errors are drawn;
// try to only use this program for linear-scale histograms,
// unless you are okay with not having error bars

R__LOAD_LIBRARY(DihBsa)
#include "Tools.h"

void RatioPlot(
  TString numerFileN="eicPlots.pythia_5x41_smear.ymin_0.03.root",
  TString denomFileN="eicPlots.pythia_5x41_smear.ymin_0.00.root",
  TString particle="hadronA",
  TString outfileN="QtRatio_PiPlus_ymin0.03",
  /*TString plotName="PperpDistLin"*/
  /*TString plotName="QtDistLin"*/
  TString plotName="QtOverQdist"
) {

  enum numerdenom {num,den};
  TFile * infile[2];
  infile[num] = new TFile(numerFileN,"READ");
  infile[den] = new TFile(denomFileN,"READ");

  const Int_t nBins = 15;
  TH1D * dist[2][nBins];
  TH1D * ratio[nBins];
  TString distN;
  TCanvas * canv[nBins];
  for(int b=0; b<nBins; b++) {
    for(int f=0; f<2; f++) {
      distN = "singlePlots/"+particle+
        "_"+plotName+"_"+Form("%d",b);
      //dist[f][b] = (TH1D*) infile[f]->Get(distN)->Clone();
      infile[f]->GetObject(distN,dist[f][b]);
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
    dist[num][b]->SetMarkerColor(kBlue);
    dist[den][b]->SetMarkerColor(kRed);
    for(int f=0; f<2; f++) {
      dist[f][b]->SetMarkerStyle(kFullCircle);
      if(plotName=="PperpDistLin")
        dist[f][b]->GetXaxis()->SetRangeUser(0,1.7);
    };
    dist[den][b]->Draw("ERR P");
    dist[num][b]->Draw("ERR P SAME");
    canv[b]->cd(2);
    canv[b]->GetPad(2)->SetGrid(1,1);
    ratio[b]->SetLineColor(kBlack);
    ratio[b]->SetMarkerColor(kBlack);
    ratio[b]->SetMarkerStyle(kFullCircle);
    ratio[b]->GetYaxis()->SetRangeUser(0,1.1);
    if(plotName=="PperpDistLin")
      ratio[b]->GetXaxis()->SetRangeUser(0,1.7);
    ratio[b]->Draw("ERR X0 P");
    canv[b]->Print(Form("%s_bin%d.pdf",outfileN.Data(),b),"pdf");
  };

  //for(int f=0; f<2; f++) infile[f]->Close();
};
