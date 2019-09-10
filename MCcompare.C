// compare MC generated and reconstructed distributions; use diagnostics.exe
// to generate plots.root files for generated set and for reconstructed set:
// - diagnostics.exe outroot.MC.gen 0x34 && mv plots{,.gen}.root
// - diagnostics.exe outroot.MC.rec 0x34 && mv plots{,.rec}.root
// - root MCcompare.C

R__LOAD_LIBRARY(DihBsa)

#include "Constants.h"
#include "EventTree.h"
#include "Tools.h"
#include "TString.h"
#include "TMath.h"

enum fenum {kGen,kRec,kRat};
int f;
TFile * infile[2];

TH1D * Fetch(Int_t f_,TString n_);
void Acceptance(TH1D ** dists_);

void MCcompare(TString fGen="plots.gen.root",
               TString fRec="plots.rec.root"
) {

  infile[kGen] = new TFile(fGen,"READ");
  infile[kRec] = new TFile(fRec,"READ");

  TH1D * WDist[3];
  TH1D * XDist[3];
  TH1D * MhDist[3];
  TH1D * ZpairDist[3];
  TH1D * xFDist[3];
  TH1D * thetaDist[3];
  TH1D * PhiHDist[3];
  TH1D * PhiRDist[3];
  TH1D * PhiHRDist[3];

  for(f=0; f<2; f++) {
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

  dists_[kRat] = (TH1D*) dists_[kRec]->Clone();

  str = "(red=gen, blue=rec) ";
  str += dists_[kGen]->GetTitle();
  dists_[kGen]->SetTitle(str);
  dists_[kRec]->SetTitle(str);
  str = "(rec/gen) ";
  str += dists_[kRat]->GetTitle();
  dists_[kRat]->SetTitle(str);

  dists_[kGen]->SetLineColor(kRed);
  dists_[kRec]->SetLineColor(kBlue);
  dists_[kRat]->SetLineColor(kBlack);

  dists_[kRat]->Divide(dists_[kGen]);
  str = "canv_"+TString(dists_[kGen]->GetName());
  TCanvas * canv = new TCanvas(str,str,1000,1000);
  canv->Divide(2,1);

  canv->cd(1);
  Int_t bigger = 
    dists_[kGen]->GetMaximum() > dists_[kRec]->GetMaximum() ? 
    kGen : kRec;
  dists_[bigger]->Draw();
  dists_[(bigger+1)%2]->Draw("same");

  canv->cd(2);
  dists_[kRat]->Draw();
};
