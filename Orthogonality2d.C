R__LOAD_LIBRARY(DihBsa)

#include "Constants.h"


void Orthogonality2d(Bool_t uniformData=0, TString infileN="ortho.root") {
  // open data hist
  TFile * infile = new TFile(infileN,"READ");
  TH2D * dataDist = (TH2D*) infile->Get("d2");

  int f,g,h,r;

  // get number of bins and ranges
  Int_t nbinsR = dataDist->GetNbinsX();
  Float_t minR = dataDist->GetXaxis()->GetXmin();
  Float_t maxR = dataDist->GetXaxis()->GetXmax();
  Float_t rangeR = maxR - minR;
  Int_t nbinsH = dataDist->GetNbinsY();
  Float_t minH = dataDist->GetYaxis()->GetXmin();
  Float_t maxH = dataDist->GetYaxis()->GetXmax();
  Float_t rangeH = maxH - minH;

  // set data weights to 1, if using uniformData
  if(uniformData) {
    for(r=1; r<=nbinsR; r++) { 
      for(h=1; h<=nbinsH; h++) { 
        dataDist->SetBinContent(r,h,6);
      };
    };
  };

  // normalization: product of ranges / integral of data hist
  Double_t normalization = rangeH*rangeR / dataDist->Integral();


  // MODULATIONS
  ////////////////////
  const Int_t NF = 8;
  TString funcT[NF];
  TString norm[NF];
  int ff=0;
  funcT[ff]="1";  norm[ff]="1/(2*PI)"; ff++;
  funcT[ff]="sin(#phi_{h}-#phi_{R})";  norm[ff]="1/(sqrt(2)*PI)"; ff++;
  funcT[ff]="sin(#phi_{R})";  norm[ff]="1/(sqrt(2)*PI)"; ff++;
  funcT[ff]="sin(#phi_{h})";  norm[ff]="1/(sqrt(2)*PI)"; ff++;
  funcT[ff]="sin(2*#phi_{h}-#phi_{R})";  norm[ff]="1/(sqrt(2)*PI)"; ff++;
  ///*
  funcT[ff]="sin(2*#phi_{h}-2*#phi_{R})";  norm[ff]="1/(sqrt(2)*PI)"; ff++;
  funcT[ff]="sin(-1*#phi_{h}+2*#phi_{R})";  norm[ff]="1/(sqrt(2)*PI)"; ff++;
  funcT[ff]="sin(3*#phi_{h}-2*#phi_{R})";  norm[ff]="1/(sqrt(2)*PI)"; ff++;
  //*/
  ////////////////////

  // define functions
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


  // |<fg>| matrix
  TH2D * orthMatrix = new TH2D("orthMatrix","|<fg>|",NF,0,NF,NF,0,NF);
  for(f=0; f<NF; f++) { 
    orthMatrix->GetXaxis()->SetBinLabel(f+1,funcT[f]);
    orthMatrix->GetYaxis()->SetBinLabel(f+1,funcT[f]);
  };


  // calculate <fg>
  TH2D * intDist[NF][NF];
  TH2D * modDist[NF];
  Float_t phiH,phiR;
  Double_t dataWeight,evalF,evalG,product,integral;
  TString intDistN,intDistT;
  TString modDistN;
  for(f=0; f<NF; f++) {
    for(g=0; g<NF; g++) {

      // generate histos, which hold the product of data*f*g
      intDistN = Form("intDist_f%d_g%d",f,g);
      intDistT = "f=" + funcT[f] + " g=" + funcT[g];
      intDist[f][g] = new TH2D(intDistN,intDistT,nbinsR,minR,maxR,nbinsH,minH,maxH);
      if(f==g) {
        modDistN = Form("modDist_f%d_g%d",f,g);
        modDist[f] = new TH2D(modDistN,funcT[f],nbinsR,minR,maxR,nbinsH,minH,maxH);
      };

      // loop over bins
      for(r=1; r<=nbinsR; r++) {
        for(h=1; h<=nbinsH; h++) {
          phiR = intDist[f][g]->GetXaxis()->GetBinCenter(r);
          phiH = intDist[f][g]->GetYaxis()->GetBinCenter(h);
          evalF = func[f]->Eval(phiR,phiH);
          evalG = func[g]->Eval(phiR,phiH);
          dataWeight = dataDist->GetBinContent(r,h);

          product = normalization * dataWeight * evalF * evalG;

          intDist[f][g]->SetBinContent(r,h,product);
          if(f==g) modDist[f]->SetBinContent(r,h,evalF);
        };
      };

      // sum data*f*g over all bins
      integral = TMath::Abs( intDist[f][g]->Integral() );
      orthMatrix->SetBinContent(f+1,g+1,integral);
    };
  };

  // get max and min for drawing
  Float_t plotMax=0;
  Float_t maxTmp;
  for(f=0; f<NF; f++) {
    for(g=0; g<NF; g++) {
      maxTmp = TMath::Max(
        TMath::Abs(intDist[f][g]->GetMaximum()),
        TMath::Abs(intDist[f][g]->GetMinimum())
      );
      plotMax = maxTmp>plotMax ? maxTmp:plotMax;
    };
  };
  for(f=0; f<NF; f++) {
    for(g=0; g<NF; g++) {
      intDist[f][g]->SetMaximum(plotMax);
      intDist[f][g]->SetMinimum(f==g?0:-plotMax);
      intDist[f][g]->SetMinimum(-plotMax);
    };
  };
      


  // draw
  TCanvas * intCanv = new TCanvas("intCanv","intCanv",1000,1000);
  TCanvas * modCanv = new TCanvas("modCanv","modCanv",1000,1000);
  Int_t pad;
  gStyle->SetOptStat(0/*1000000*/);
  intCanv->Divide(NF,NF);
  modCanv->Divide(NF,NF);
  for(f=0; f<NF; f++) {
    for(g=0; g<NF; g++) {
      pad = (NF-1-g)*NF+f+1;
      intCanv->cd(pad);
      intDist[f][g]->Draw("colz");
      if(f==g) {
        modCanv->cd(pad);
        modDist[f]->Draw("colz");
      };
    };
  };

  TCanvas * matCanv = new TCanvas("matCanv","matCanv",1000,1000);
  orthMatrix->Draw("colztext");

};
