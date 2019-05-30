R__LOAD_LIBRARY(DihBsa)

#include "Constants.h"

void LookForPi0s(TString dir="outroot") {
  TString files = dir + "/*.root";
  TChain * c = new TChain("tree");
  c->Add(files);

  TString cutStr;

  cutStr = Form("pairType==%d",pairP0);
  cutStr += " && diphE>2";
  cutStr += " && diphZ<0.6";

  printf("cutStr = %s\n",cutStr.Data());

  new TCanvas();
  c->Draw("diphM>>hMass(200,0,1)",cutStr);

  new TCanvas();
  c->Draw("diphAlpha>>hAlpha(200,0,2)",cutStr);

  new TCanvas();
  c->Draw("diphAlpha:diphM>>hAlphaM(200,0,1,200,0,2)",cutStr,"colz");

  new TCanvas();
  c->Draw("diphM:diphE>>hME(200,0,5,200,0,1)",cutStr,"colz");

  new TCanvas();
  c->Draw("diphM:diphZ>>hMZ(200,0,1,200,0,1)",cutStr,"colz");

  new TCanvas();
  c->Draw("diphM:diphPt>>hMPt(200,0,1,200,0,1)",cutStr,"colz");
};

  

