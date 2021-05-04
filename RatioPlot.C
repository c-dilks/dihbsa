// plot ratio of specific plots in eicPlots*.root files
//
// WARNING: drawing histograms with variable bin widths (e.g.,
// even bins in log scale), does not work if errors are drawn;
// try to only use this program for linear-scale histograms,
// unless you are okay with not having error bars

R__LOAD_LIBRARY(DihBsa)
#include "Tools.h"

void RatioPlot(
  TString numerFileN="roots_pTgt1.0/ymin_0.05.root",
  TString denomFileN="roots_pTgt1.0/ymin_0.00.root",
  TString particle="hadronA",
  TString outfileN="Pperp_PiPlus_ymin0.05",
  TString plotName="PperpDistLin"
  /*TString plotName="QtDistLin"*/
  /*TString plotName="QtOverQdist"*/
) {

  enum numerdenom {num,den};
  TFile * infile[2];
  infile[num] = new TFile(numerFileN,"READ");
  infile[den] = new TFile(denomFileN,"READ");

  Bool_t plotRatioOnly = 0;

  const Int_t nBins = 15;
  TH1D * dist[2][nBins];
  TH1D * ratio[nBins];
  Double_t norm;
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
      dist[f][b]->Rebin(4);
      norm = dist[f][b]->Integral();
      //dist[f][b]->Scale(1./norm);
    };

    ratio[b] = (TH1D*) dist[num][b]->Clone();
    ratio[b]->Divide(dist[den][b]);

    if(plotRatioOnly)
      canv[b] = new TCanvas(Form("canv%d",b),Form("canv%d",b), 800,700);
    else
      canv[b] = new TCanvas(Form("canv%d",b),Form("canv%d",b), 1600,700);

    if(!plotRatioOnly) {
      canv[b]->Divide(2,1);
      canv[b]->cd(1);
      canv[b]->GetPad(1)->SetGrid(1,1);
      //canv[b]->GetPad(1)->SetLogx();
      //canv[b]->GetPad(1)->SetLogy();
      canv[b]->GetPad(1)->SetBottomMargin(0.15);
      canv[b]->GetPad(1)->SetLeftMargin(0.15);
      dist[num][b]->SetLineColor(kBlue);
      dist[den][b]->SetLineColor(kRed);
      dist[num][b]->SetMarkerColor(kBlue);
      dist[den][b]->SetMarkerColor(kRed);
      for(int f=0; f<2; f++) {
        dist[f][b]->SetMarkerStyle(kFullCircle);
        dist[f][b]->SetMarkerSize(1.5);
        dist[f][b]->SetLineWidth(2);
        dist[f][b]->GetXaxis()->SetLabelSize(0.06);
        dist[f][b]->GetYaxis()->SetLabelSize(0.06);
        dist[f][b]->GetXaxis()->SetTitleSize(0.06);
        dist[f][b]->GetYaxis()->SetTitleSize(0.06);
        dist[f][b]->GetXaxis()->SetTitleOffset(1.2);
        if(plotName=="PperpDistLin")
          dist[f][b]->GetXaxis()->SetRangeUser(0,1.7);
      };
      dist[den][b]->Draw("ERR P");
      dist[num][b]->Draw("ERR P SAME");
      canv[b]->cd(2);
      canv[b]->GetPad(2)->SetGrid(1,1);
      canv[b]->GetPad(2)->SetBottomMargin(0.15);
      canv[b]->GetPad(2)->SetLeftMargin(0.15);
    } else {
      canv[b]->cd();
      canv[b]->SetGrid(1,1);
      canv[b]->SetBottomMargin(0.15);
      canv[b]->SetLeftMargin(0.15);
    };
    ratio[b]->GetYaxis()->SetTitle("ratio( y>ymin / total )");
    ratio[b]->SetLineColor(kBlack);
    ratio[b]->SetMarkerColor(kBlack);
    ratio[b]->SetMarkerStyle(kFullCircle);
    ratio[b]->SetMarkerSize(1.5);
    ratio[b]->SetLineWidth(2);
    ratio[b]->GetXaxis()->SetLabelSize(0.06);
    ratio[b]->GetYaxis()->SetLabelSize(0.06);
    ratio[b]->GetXaxis()->SetTitleSize(0.06);
    ratio[b]->GetYaxis()->SetTitleSize(0.06);
    ratio[b]->GetXaxis()->SetTitleOffset(1.2);
    ratio[b]->GetYaxis()->SetRangeUser(0,1.3);
    if(plotName=="PperpDistLin") ratio[b]->GetXaxis()->SetRangeUser(0,1.7);
    ratio[b]->Draw("ERR X0 P");
    canv[b]->Print(Form("%s_bin%d.png",outfileN.Data(),b),"png");
  };

  //for(int f=0; f<2; f++) infile[f]->Close();
};
