// compare different fiducial cuts; use diagnostics.exe
// to generate required plots.*.root files that you want to compare

R__LOAD_LIBRARY(DihBsa)

#include "Constants.h"
#include "EventTree.h"
#include "Tools.h"
#include "TString.h"
#include "TMath.h"

enum fenum {kCut0,kCut1,kCut2,N};
enum renum {kCut0Cut1=N,kCut0Cut2,NR};
int f;
TFile * infile[N];
Int_t plotColor[NR];
TString name[N];

TH1D * Fetch(Int_t f_,TString n_);
void CompareDists(TH1D ** dists_);

void FiducialCompare(TString fCut0="plots.fidNone.root",
  TString fCut1="plots.fidLoosePCALandDC.root",TString fCut2="plots.fidTightPCALandDC.root"
  /*
  TString fCut1="plots.fidLoosePCAL.root",TString fCut2="plots.fidTightPCAL.root"
  TString fCut1="plots.fidLooseDC.root",TString fCut2="plots.fidTightDC.root"
  */
) {

  name[kCut0] = fCut0;
  name[kCut1] = fCut1;
  name[kCut2] = fCut2;
  for(f=0; f<N; f++) {
    name[f](TRegexp("^.*plots\\."))="";
    name[f](TRegexp("\\.root$"))="";
  };

  infile[kCut0] = new TFile(fCut0,"READ");
  infile[kCut1] = new TFile(fCut1,"READ");
  infile[kCut2] = new TFile(fCut2,"READ");

  TH1D * Q2Dist[NR];
  //TH1D * WDist[NR]; // WDist in diagnostics.cpp does not have fiducial cuts applied
  TH1D * XDist[NR];
  TH1D * MhDist[NR];
  TH1D * ZpairDist[NR];
  TH1D * xFDist[NR];
  TH1D * thetaDist[NR];
  TH1D * PhiHDist[NR];
  TH1D * PhiRDist[NR];
  TH1D * PhiHRDist[NR];

  for(f=0; f<N; f++) {
    Q2Dist[f] = Fetch(f,"Q2Dist");
    //WDist[f] = Fetch(f,"WDist");
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

  plotColor[kCut0] = kBlack;
  plotColor[kCut1] = kOrange+1;
  plotColor[kCut2] = kBlue;
  plotColor[kCut0Cut1] = kBlack;
  plotColor[kCut0Cut2] = kBlack;

  CompareDists(Q2Dist);
  //CompareDists(WDist);
  CompareDists(XDist);
  CompareDists(MhDist);
  CompareDists(ZpairDist);
  CompareDists(xFDist);
  CompareDists(thetaDist);
  CompareDists(PhiHDist);
  CompareDists(PhiRDist);
  CompareDists(PhiHRDist);

};

  
TH1D * Fetch(Int_t f_,TString n_) {
  TH1D * d_ = (TH1D*) infile[f_]->Get(n_);
  Double_t ent = d_->GetEntries();
  d_->Scale(1/ent);
  return d_;
};


void CompareDists(TH1D ** dists_) {
  TString str;

  dists_[kCut0Cut1] = (TH1D*) dists_[kCut0]->Clone();
  dists_[kCut0Cut2] = (TH1D*) dists_[kCut0]->Clone();

  str = "[ "+name[kCut0]+" / "+name[kCut1]+" ] -- " + 
        TString(dists_[kCut0Cut1]->GetTitle());
  dists_[kCut0Cut1]->SetTitle(str);
  str = "[ "+name[kCut0]+" / "+name[kCut2]+" ] -- " + 
        TString(dists_[kCut0Cut2]->GetTitle());
  dists_[kCut0Cut2]->SetTitle(str);

  for(f=0; f<NR; f++) {
    dists_[f]->SetLineWidth(2);
    dists_[f]->SetMarkerSize(0.7);
    dists_[f]->SetLineColor(plotColor[f]);
    dists_[f]->SetMarkerColor(plotColor[f]);
    dists_[f]->SetMarkerStyle(kFullCircle);
  };
  dists_[kCut0]->SetMarkerStyle(kPlus);

  TLegend * leg[NR];
  for(f=0; f<NR; f++) leg[f] = new TLegend(0.6, 0.8, 0.9, 0.9);
  leg[kCut0Cut1]->AddEntry(dists_[kCut0],name[kCut0],"PLE");
  leg[kCut0Cut1]->AddEntry(dists_[kCut1],name[kCut1],"PLE");
  leg[kCut0Cut2]->AddEntry(dists_[kCut0],name[kCut0],"PLE");
  leg[kCut0Cut2]->AddEntry(dists_[kCut2],name[kCut2],"PLE");


  dists_[kCut0Cut1]->Divide(dists_[kCut1]);
  dists_[kCut0Cut2]->Divide(dists_[kCut2]);

  str = "canv_"+TString(dists_[kCut1]->GetName());
  TCanvas * canv = new TCanvas(str,str,1000,1000);
  canv->Divide(2,2);
  Int_t d[2];
  Float_t minFit, maxFit;

  canv->cd(1);
  if(dists_[kCut0]->GetMaximum() > dists_[kCut1]->GetMaximum()) { d[0]=kCut0; d[1]=kCut1; }
  else { d[0]=kCut1; d[1]=kCut0; };
  dists_[d[0]]->Draw();
  dists_[d[1]]->Draw("same");
  dists_[kCut0]->Draw("same");
  leg[kCut0Cut1]->Draw();

  canv->cd(2);
  dists_[kCut0Cut1]->Draw();
  minFit = dists_[kCut0Cut1]->GetXaxis()->GetXmin();
  maxFit = dists_[kCut0Cut1]->GetXaxis()->GetXmax();
  dists_[kCut0Cut1]->Fit("pol0","","",minFit,maxFit);

  canv->cd(3);
  if(dists_[kCut0]->GetMaximum() > dists_[kCut2]->GetMaximum()) { d[0]=kCut0; d[1]=kCut2; }
  else { d[0]=kCut2; d[1]=kCut0; };
  dists_[d[0]]->Draw();
  dists_[d[1]]->Draw("same");
  dists_[kCut0]->Draw("same");
  leg[kCut0Cut2]->Draw();

  canv->cd(4);
  dists_[kCut0Cut2]->Draw();
  minFit = dists_[kCut0Cut2]->GetXaxis()->GetXmin();
  maxFit = dists_[kCut0Cut2]->GetXaxis()->GetXmax();
  dists_[kCut0Cut2]->Fit("pol0","","",minFit,maxFit);
};
