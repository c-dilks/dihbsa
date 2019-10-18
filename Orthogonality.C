R__LOAD_LIBRARY(DihBsa)

#include "Constants.h"


void Orthogonality(TString infileN="plots.data.root") {
  TFile * infile = new TFile(infileN,"READ");
  TH2D * dataDist = (TH2D*) infile->Get("PhiHvsPhiR");

  int f,g,h,r;

  Int_t nbinsR = dataDist->GetNbinsX();
  Float_t minR = dataDist->GetXaxis()->GetXmin();
  Float_t maxR = dataDist->GetXaxis()->GetXmax();
  Int_t nbinsH = dataDist->GetNbinsY();
  Float_t minH = dataDist->GetYaxis()->GetXmin();
  Float_t maxH = dataDist->GetYaxis()->GetXmax();

  /*
  // TEST: set data weights to 1
  for(r=1; r<=nbinsR; r++) { for(h=1; h<=nbinsH; h++) { 
    dataDist->SetBinContent(r,h,1);
  }};
  */

  dataDist->Scale( (nbinsR*nbinsH) / dataDist->Integral() );

  const Int_t NF = 5;
  TString funcT[NF];
  TString norm[NF];
  int ff=0;
  funcT[ff]="1";  norm[ff]="1/(2*PI)"; ff++;
  funcT[ff]="sin(#phi_{h}-#phi_{R})";  norm[ff]="1/(sqrt(2)*PI)"; ff++;
  funcT[ff]="sin(#phi_{R})";  norm[ff]="1/(sqrt(2)*PI)"; ff++;
  funcT[ff]="sin(#phi_{h})";  norm[ff]="1/(sqrt(2)*PI)"; ff++;
  funcT[ff]="sin(#phi_{h}+#phi_{R})";  norm[ff]="1/(sqrt(2)*PI)"; ff++;

  TString formu[NF];
  TString funcN[NF];
  TF2 * func[NF];
  printf("\nFUNCTIONS:\n");
  for(f=0; f<NF; f++) {
    formu[f] = funcT[f];
    formu[f](TRegexp("sin")) = "TMath::Sin";
    formu[f](TRegexp("#phi_{h}")) = "y";
    formu[f](TRegexp("#phi_{R}")) = "x";

    norm[f](TRegexp("PI")) = "TMath::Pi()";
    norm[f](TRegexp("sqrt")) = "TMath::Sqrt";
    formu[f] = "("+norm[f]+")*("+formu[f]+")";

    printf("%s -> %s\n",funcT[f].Data(),formu[f].Data());
    funcN[f] = Form("func%d",f);
    func[f] = new TF2(funcN[f],formu[f],minR,maxR,minH,maxH);
  };

  TF2 * productFunc[NF][NF];
  TString productFuncN;
  Double_t anInt[NF][NF];
  for(f=0; f<NF; f++) {
    for(g=0; g<NF; g++) {
      productFuncN = Form("productFunc_f%d_g%d",f,g);
      productFunc[f][g] = new TF2(productFuncN,
        TString("("+formu[f]+")*("+formu[g]+")"),minR,maxR,minH,maxH);
      anInt[f][g] = productFunc[f][g]->Integral(minR,maxR,minH,maxH);
    };
  };


  TH2D * orthMatrix = new TH2D("orthMatrix","|<fg>|",NF,0,NF,NF,0,NF);

  TH2D * intDist[NF][NF];
  Float_t phiH,phiR;
  Double_t dataWeight,evalF,evalG,product,integral;
  TString intDistN,intDistT;
  for(f=0; f<NF; f++) {
    for(g=0; g<NF; g++) {
      intDistN = Form("intDist_f%d_g%d",f,g);
      intDistT = "f=" + funcT[f] + " g=" + funcT[g];
      intDist[f][g] = new TH2D(intDistN,intDistT,nbinsR,minR,maxR,nbinsH,minH,maxH);
      for(r=1; r<=nbinsR; r++) {
        for(h=1; h<=nbinsH; h++) {
          phiR = intDist[f][g]->GetXaxis()->GetBinCenter(r);
          phiH = intDist[f][g]->GetYaxis()->GetBinCenter(h);
          evalF = func[f]->Eval(phiR,phiH);
          evalG = func[g]->Eval(phiR,phiH);
          dataWeight = dataDist->GetBinContent(r,h);

          product = dataWeight * evalF * evalG;

          intDist[f][g]->SetBinContent(r,h,product);
        };
      };
      integral = TMath::Abs( intDist[f][g]->Integral("width") );
      orthMatrix->SetBinContent(f+1,g+1,integral);
    };
  };

  TCanvas * intCanv = new TCanvas("intCanv","intCanv",1000,1000);
  gStyle->SetOptStat(0/*1000000*/);
  intCanv->Divide(NF,NF);
  for(f=0; f<NF; f++) {
    for(g=0; g<NF; g++) {
      intCanv->cd((NF-1-g)*NF+f+1);
      intDist[f][g]->Draw("colz");
    };
  };

  TCanvas * matCanv = new TCanvas("matCanv","matCanv",1000,1000);
  orthMatrix->Draw("colztext");

};
