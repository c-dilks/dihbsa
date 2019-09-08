R__LOAD_LIBRARY(DihBsa)

#include "TString.h"
#include "TMultiGraph.h"
#include "TRegexp.h"

void DrawMultiGr(TString infileN="spinFinal.root",TString suffix="M") {
  TFile * infile = new TFile(infileN,"READ");

  TString runs = infileN;
  runs(TRegexp(".root$")) = "";
  runs(TRegexp("^spinFinal")) = "";
  runs(TRegexp("^.")) = "";

  TString pngname = "multigr."+runs+".png";
  printf("pngname = %s\n",pngname.Data());

  TCanvas * c = new TCanvas("c","c",600,600);
  c->SetGrid(1,1);
  TMultiGraph * gr = (TMultiGraph*) infile->Get(TString("multiGr_"+suffix));
  TString grT = gr->GetTitle();
  grT = "runs " + runs + " -- " + grT;
  gr->SetTitle(grT);

  gr->Draw("LAPE");
  gr->GetYaxis()->SetRangeUser(-0.03,0.05);
  gr->Draw("LAPE");


  c->Print(pngname,"png");
};


