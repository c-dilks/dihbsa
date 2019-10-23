R__LOAD_LIBRARY(DihBsa)

#include "Constants.h"

Bool_t disableLegendre;
Double_t Azimod(Int_t twist_, Int_t m_, Float_t phiH_, Float_t phiR_);
Double_t Legendre(Int_t l_, Int_t m_, Float_t theta_);

// weightSetting: 
// 0 = use acceptance from infile (ortho.root)
// 1 = uniform across all variables
// 2+ = see switch statement

void Orthogonality3d(Int_t weightSetting=0, TString infileN="ortho.root") {

  
  // OPTIONS
  ///////////////////
  disableLegendre = 1;
  Int_t LMAX = 2;
  ///////////////////


  // open data hist
  TFile * infile = new TFile(infileN,"READ");
  TH3D * dataDist = (TH3D*) infile->Get("d3");

  int f,g,h,r,t;

  // get number of bins and bin widths
  Int_t nbinsR = dataDist->GetNbinsX();
  Float_t minR = dataDist->GetXaxis()->GetXmin();
  Float_t maxR = dataDist->GetXaxis()->GetXmax();
  //Float_t rangeR = maxR - minR;
  //Float_t widthR = dataDist->GetXaxis()->GetBinWidth(0);

  Int_t nbinsH = dataDist->GetNbinsY();
  Float_t minH = dataDist->GetYaxis()->GetXmin();
  Float_t maxH = dataDist->GetYaxis()->GetXmax();
  //Float_t rangeH = maxH - minH;
  //Float_t widthH = dataDist->GetYaxis()->GetBinWidth(0);

  Int_t nbinsT = dataDist->GetNbinsZ();
  Float_t minT = dataDist->GetZaxis()->GetXmin();
  Float_t maxT = dataDist->GetZaxis()->GetXmax();
  //Float_t rangeT = maxT - minT;
  //Float_t widthT = dataDist->GetZaxis()->GetBinWidth(0);

  //Double_t binSize = widthR * widthH * widthT;


  // set test weight distributions (if weightSetting==0, just uses those from infile) 
  if(weightSetting>0) {
    for(r=1; r<=nbinsR; r++) { 
      for(h=1; h<=nbinsH; h++) { 
        for(t=1; t<=nbinsT; t++) { 
          switch(weightSetting) {
            case 1: dataDist->SetBinContent(r,h,t,
              1
            ); break;
            case 2: dataDist->SetBinContent(r,h,t,
              r < (Double_t)nbinsR/2. ? 1:0 
            ); break;
            case 3: dataDist->SetBinContent(r,h,t,
              h < (Double_t)nbinsH/2. ? 1:0 
            ); break;
            case 4: dataDist->SetBinContent(r,h,t,
              t < (Double_t)nbinsT/2. ? 1:0 
            ); break;
            default: fprintf(stderr,"ERROR: bad weightSetting\n"); return;
          };
        };
      };
    };
  };


  // generate list of (l,m,twist) values
  Int_t NFi = 0;
  enum idx_enum{kL,kM,kTwist,Nidx};
  Int_t idx[100][Nidx];
  Int_t l,m,twist;
  idx[NFi][kL]=-1; idx[NFi][kM]=0; idx[NFi][kTwist]=0; NFi++; // (constant)
  for(l=0; l<=LMAX; l++) {
    for(m=0; m<=l; m++) {
      for(twist=2; twist<=3; twist++) {
        if(disableLegendre && l<LMAX) continue;
        if((twist==2 && m>0) || twist==3) {
          idx[NFi][kL]=l; idx[NFi][kM]=m; idx[NFi][kTwist]=twist; NFi++;
          if(twist==3 && m>0) {
            idx[NFi][kL]=l; idx[NFi][kM]=-m; idx[NFi][kTwist]=twist; NFi++;
          };
        } 
      };
    };
  };
  const Int_t NF = NFi;
        


  // |<fg>| matrix
  TH2D * orthMatrix = new TH2D("orthMatrix","|<fg>|",NF,0,NF,NF,0,NF);
  TString funcT[NF];
  for(f=0; f<NF; f++) { 
    if(!disableLegendre) 
      funcT[f] = Form("|%d,%d>_{%d}",idx[f][kL],idx[f][kM],idx[f][kTwist]);
    else
      funcT[f] = Form("|L,%d>_{%d}",idx[f][kM],idx[f][kTwist]);
    if(idx[f][kL]==-1) funcT[f]="const";
    printf("%s\n",funcT[f].Data());
    orthMatrix->GetXaxis()->SetBinLabel(f+1,funcT[f]);
    orthMatrix->GetYaxis()->SetBinLabel(f+1,funcT[f]);
  };


  // calculate <fg>
  TH3D * intDist[NF][NF];
  Float_t phiH,phiR,theta;
  Double_t dataWeight,evalF,evalG,product,integral;
  TString intDistN,intDistT;
  Double_t weightedNorm[NF];
  for(f=0; f<NF; f++) {
    for(g=0; g<NF; g++) {

      // generate histos, which hold the product of data*f*g
      intDistN = Form("intDist_f%d_g%d",f,g);
      intDistT = "f: " + funcT[f] + ",  g: " + funcT[g];
      intDist[f][g] = new TH3D(intDistN,intDistT,
        nbinsR,minR,maxR,
        nbinsH,minH,maxH,
        nbinsT,minT,maxT);

      // loop over bins
      for(r=1; r<=nbinsR; r++) {
        for(h=1; h<=nbinsH; h++) {
          for(t=1; t<=nbinsT; t++) {
            phiR = intDist[f][g]->GetXaxis()->GetBinCenter(r);
            phiH = intDist[f][g]->GetYaxis()->GetBinCenter(h);
            theta = intDist[f][g]->GetZaxis()->GetBinCenter(t);

            evalF = Legendre( idx[f][kL], idx[f][kM], theta);
            evalG = Legendre( idx[g][kL], idx[g][kM], theta);

            evalF *= Azimod( idx[f][kTwist], idx[f][kM], phiH, phiR);
            evalG *= Azimod( idx[g][kTwist], idx[g][kM], phiH, phiR);

            dataWeight = dataDist->GetBinContent(r,h,t);

            product = dataWeight * evalF * evalG;

            intDist[f][g]->SetBinContent(r,h,t,product);
          };
        };
      };

      // computed weighted normalization for function f as 1/sqrt(<ff>)
      if(f==g) {
        integral = TMath::Abs( intDist[f][g]->Integral("width") );
        weightedNorm[f] = 1 / TMath::Sqrt(integral);
      };
    };
  };



  // fill orthogonality matrix
  for(f=0; f<NF; f++) {
    for(g=0; g<NF; g++) {

      // normalize data*f*g by 1/sqrt(<ff><gg>)
      intDist[f][g]->Scale(weightedNorm[f]*weightedNorm[g]);

      // sum data*f*g over all bins
      integral = TMath::Abs( intDist[f][g]->Integral("width") );
      orthMatrix->SetBinContent(f+1,g+1,integral);
    };
  };

  // draw
  gStyle->SetOptStat(0);
  gStyle->SetPalette(kBird);
  gStyle->SetPaintTextFormat(".3f");
  TCanvas * matCanv = new TCanvas("matCanv","matCanv",1000,1000);
  orthMatrix->Draw("colztext");

};


Double_t Azimod(Int_t twist_, Int_t m_, Float_t phiH_, Float_t phiR_) {

  Double_t retval_;

  /*
  Double_t norm_ = 1 / ( TMath::Sqrt(2) * PI ); //+++
  if(twist_==0) return 1 / ( 2*PI ); // (constant) //+++
  */
  Double_t norm_ = 1; // normalization not needed

  if(twist_==0) return 1; // (constant)
  else if(twist_==2) retval_ = TMath::Sin( m_ * (phiH_-phiR_) );
  else if(twist_==3) retval_ = TMath::Sin( (1-m_)*phiH_ + m_*phiR_ );
  else {
    fprintf(stderr,"ERROR: bad twist\n");
    return 0;
  };

  return norm_*retval_;
};


Double_t Legendre(Int_t l_, Int_t m_, Float_t theta_) {
  Double_t retval_, norm_;
  m_ = TMath::Abs(m_);


  if(l_<0 || disableLegendre) return 1 / TMath::Sqrt(PI); // (constant)
  if(m_<0||m_>l_) {
    fprintf(stderr,"ERROR: bad m\n");
    return 0;
  }
  if(l_==0) {
    switch(m_) {
      case 0:
        retval_ = 1;
        norm_ = 1 / TMath::Sqrt(PI);
        break;
    };
  } else if(l_==1) {
    switch(m_) {
      case 0:
        retval_ = TMath::Cos(theta_);
        norm_ = TMath::Sqrt( 2/PI );
        break;
      case 1:
        retval_ = TMath::Sin(theta_);
        norm_ = TMath::Sqrt( 2/PI );
        break;
    };
  } else if(l_==2) {
    switch(m_) {
      case 0:
        retval_ = 0.5*( 3*TMath::Power(TMath::Cos(theta_),2) - 1 );
        norm_ = 4 * TMath::Sqrt( 2/(11*PI) );
        break;
      case 1:
        retval_ = 3 * TMath::Sin(theta_) * TMath::Cos(theta_);
        norm_ = (2./3.) * TMath::Sqrt( 2/PI );
        break;
      case 2:
        retval_ = 3 * TMath::Power(TMath::Sin(theta_),2);
        norm_ = (2./3.) * TMath::Sqrt( 2/(3*PI) );
        break;
    };
  } else {
    fprintf(stderr,"ERROR: bad l\n");
    return 0;
  };

  norm_ = 1; // normalization not needed
  
  return retval_ * norm_;
};
