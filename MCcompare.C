// compare MC generated and reconstructed distributions; use diagnostics.exe
// to generate plots.root files for generated, reconstructed, and data sets:
// - diagnostics.exe outroot.MC.gen 0x34 && mv plots{,.gen}.root
// - diagnostics.exe outroot.MC.rec 0x34 && mv plots{,.rec}.root
// - diagnostics.exe outroot.fall18.some 0x34 && mv plots{,.data}.root
// - root MCcompare.C

R__LOAD_LIBRARY(DihBsa)

#include "Constants.h"
#include "EventTree.h"
#include "Tools.h"
#include "TString.h"
#include "TMath.h"

enum fenum {kGen,kRec,kData,N};
enum renum {kRecGen=N,kRecData,NR};
int f;
TFile * infile[N];
Int_t plotColor[NR];

TH1D * Fetch(Int_t f_,TString n_);
void Acceptance(TH1D ** dists_);

void MCcompare(TString fGen="plots.gen.root",
               TString fRec="plots.rec.root",
               TString fData="plots.data.root"
) {

  infile[kGen] = new TFile(fGen,"READ");
  infile[kRec] = new TFile(fRec,"READ");
  infile[kData] = new TFile(fData,"READ");

  TH1D * WDist[NR];
  TH1D * XDist[NR];
  TH1D * MhDist[NR];
  TH1D * ZpairDist[NR];
  TH1D * xFDist[NR];
  TH1D * thetaDist[NR];
  TH1D * PhiHDist[NR];
  TH1D * PhiRDist[NR];
  TH1D * PhiHRDist[NR];

  for(f=0; f<N; f++) {
    WDist[f] = Fetch(f,"WDist");
    XDist[f] = Fetch(f,"XDist");
    MhDist[f] = Fetch(f,"MhDist");
    ZpairDist[f] = Fetch(f,"ZpairDist");
    xFDist[f] = Fetch(f,"xFDist");
    thetaDist[f] = Fetch(f,"thetaDist");
    PhiHDist[f] = Fetch(f,"PhiHDist");
    PhiRDist[f] = Fetch(f,"PhiRDist");
    PhiHRDist[f] = Fetch(f,"PhiHRDist");
  };

  gStyle->SetOptStat(0);

  plotColor[kGen] = kRed;
  plotColor[kRec] = kBlue;
  plotColor[kData] = kMagenta;
  plotColor[kRecGen] = kBlack;
  plotColor[kRecData] = kBlack;

  Acceptance(WDist);
  Acceptance(XDist);
  Acceptance(MhDist);
  Acceptance(ZpairDist);
  Acceptance(xFDist);
  Acceptance(thetaDist);
  Acceptance(PhiHDist);
  Acceptance(PhiRDist);
  Acceptance(PhiHRDist);

};

  
TH1D * Fetch(Int_t f_,TString n_) {
  TH1D * d_ = (TH1D*) infile[f_]->Get(n_);
  Double_t ent = d_->GetEntries();
  d_->Scale(1/ent);
  return d_;
};


void Acceptance(TH1D ** dists_) {
  TString str;

  dists_[kRecGen] = (TH1D*) dists_[kRec]->Clone();
  dists_[kRecData] = (TH1D*) dists_[kRec]->Clone();

  /*
  str = "(red=gen, blue=rec, green=data) ";
  str += dists_[kGen]->GetTitle();
  for(f=0; f<N; f++) dists_[f]->SetTitle(str);
  */

  str = "[ MCrec / MCgen ] -- " + TString(dists_[kRecGen]->GetTitle());
  dists_[kRecGen]->SetTitle(str);
  str = "[ MCrec / data ] -- " + TString(dists_[kRecData]->GetTitle());
  dists_[kRecData]->SetTitle(str);

  for(f=0; f<NR; f++) {
    dists_[f]->SetLineWidth(2);
    dists_[f]->SetMarkerSize(0.7);
    dists_[f]->SetLineColor(plotColor[f]);
    dists_[f]->SetMarkerColor(plotColor[f]);
    dists_[f]->SetMarkerStyle(f==kRec?kOpenCircle:kFullCircle);
  };

  TLegend * leg[NR];
  for(f=0; f<NR; f++) leg[f] = new TLegend(0.6, 0.8, 0.9, 0.9);
  leg[kRecGen]->AddEntry(dists_[kRec],"MC reconstructed","PLE");
  leg[kRecGen]->AddEntry(dists_[kGen],"MC generated","PLE");
  leg[kRecData]->AddEntry(dists_[kRec],"MC reconstructed","PLE");
  leg[kRecData]->AddEntry(dists_[kData],"data","PLE");


  dists_[kRecGen]->Divide(dists_[kGen]);
  dists_[kRecData]->Divide(dists_[kData]);

  str = "canv_"+TString(dists_[kGen]->GetName());
  //TCanvas * canv = new TCanvas(str,str,1000,1000);
  TCanvas * canv = new TCanvas(str,str,1200,400);
  //canv->Divide(2,2);
  canv->Divide(3,1);
  Int_t d[2];

  canv->cd(1);
  if(dists_[kRec]->GetMaximum() > dists_[kGen]->GetMaximum()) { d[0]=kRec; d[1]=kGen; }
  else { d[0]=kGen; d[1]=kRec; };
  dists_[d[0]]->Draw();
  dists_[d[1]]->Draw("same");
  leg[kRecGen]->Draw();

  //canv->cd(2);
  //dists_[kRecGen]->Draw();

  //canv->cd(3);
  canv->cd(2);
  if(dists_[kRec]->GetMaximum() > dists_[kData]->GetMaximum()) { d[0]=kRec; d[1]=kData; }
  else { d[0]=kData; d[1]=kRec; };
  dists_[d[0]]->Draw();
  dists_[d[1]]->Draw("same");
  leg[kRecData]->Draw();

  //canv->cd(4);
  canv->cd(3);
  dists_[kRecData]->Draw();
  Float_t minFit = dists_[kRecData]->GetXaxis()->GetXmin();
  Float_t maxFit = dists_[kRecData]->GetXaxis()->GetXmax();
  dists_[kRecData]->Fit("pol0","","",minFit,maxFit);
};
