R__LOAD_LIBRARY(DihBsa)

#include "Constants.h"

Double_t Azimod(Int_t twist_, Int_t m_, Float_t phiH_, Float_t phiR_);
Double_t Legendre(Int_t l_, Int_t m_, Float_t theta_);

void Orthogonality3d(Bool_t uniformData=0, TString infileN="ortho.root") {
  // open data hist
  TFile * infile = new TFile(infileN,"READ");
  TH3D * dataDist = (TH3D*) infile->Get("d3");

  int f,g,h,r,t;

  // get number of bins and bin widths
  Int_t nbinsR = dataDist->GetNbinsX();
  Float_t minR = dataDist->GetXaxis()->GetXmin();
  Float_t maxR = dataDist->GetXaxis()->GetXmax();
  //Float_t rangeR = maxR - minR;
  Float_t widthR = dataDist->GetXaxis()->GetBinWidth(0);

  Int_t nbinsH = dataDist->GetNbinsY();
  Float_t minH = dataDist->GetYaxis()->GetXmin();
  Float_t maxH = dataDist->GetYaxis()->GetXmax();
  //Float_t rangeH = maxH - minH;
  Float_t widthH = dataDist->GetYaxis()->GetBinWidth(0);

  Int_t nbinsT = dataDist->GetNbinsZ();
  Float_t minT = dataDist->GetZaxis()->GetXmin();
  Float_t maxT = dataDist->GetZaxis()->GetXmax();
  //Float_t rangeT = maxT - minT;
  Float_t widthT = dataDist->GetZaxis()->GetBinWidth(0);

  Double_t binSize = widthR * widthH * widthT;


  // set data weights to 1, if using uniformData
  if(uniformData) {
    for(r=1; r<=nbinsR; r++) { 
      for(h=1; h<=nbinsH; h++) { 
        for(t=1; t<=nbinsT; t++) { 
          //if( r < (Double_t)nbinsR/2. /*&& h < (Double_t)nbinsH/2.*/)
            dataDist->SetBinContent(r,h,t,1);
          //else
            //dataDist->SetBinContent(r,h,t,0);
        };
      };
    };
  };

  // data normalization: product of ranges / integral of data hist
  Int_t nFilledBins = 0;
  for(r=1; r<=nbinsR; r++) { 
    for(h=1; h<=nbinsH; h++) { 
      for(t=1; t<=nbinsT; t++) { 
        if(dataDist->GetBinContent(r,h,t)>0) nFilledBins++;
      };
    };
  };
  //Double_t dataNorm = nFilledBins / dataDist->Integral();
  //Double_t dataNorm = nbinsR*nbinsH*nbinsT / dataDist->Integral();
  //Double_t dataNorm = 1 / dataDist->Integral();
  //printf("dataDist %.2f%% filled\n",100.*nFilledBins/(nbinsR*nbinsH*nbinsT));
  //dataDist->Scale(dataNorm);


  // generate list of (l,m,twist) values
  Int_t NFi = 0;
  enum idx_enum{kL,kM,kTwist,Nidx};
  Int_t idx[100][Nidx];
  Int_t l,m,twist;
  Int_t LMAX = 1;
  idx[NFi][kL]=-1; idx[NFi][kM]=0; idx[NFi][kTwist]=0; NFi++; // (constant)
  for(l=0; l<=LMAX; l++) {
    for(twist=2; twist<=3; twist++) {
      for(m=0; m<=l; m++) {
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
    funcT[f] = Form("|%d,%d>_{%d}",idx[f][kL],idx[f][kM],idx[f][kTwist]);
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
      nFilledBins = 0;
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

            if(TMath::Abs(product)>0.0000000001) nFilledBins++;
            //printf("p=%f\n",product);

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

      //intDist[f][g]->Scale(nFilledBins*dataNorm);
      //printf("    %.2f%% filled\n",100.*nFilledBins/(nbinsR*nbinsH*nbinsT));

      // normalize data*f*g by weighted normalizations 1/sqrt(<ff><gg>)
      intDist[f][g]->Scale(weightedNorm[f]*weightedNorm[g]);

      // sum data*f*g over all bins
      integral = TMath::Abs( intDist[f][g]->Integral("width") );
      orthMatrix->SetBinContent(f+1,g+1,integral);
    };
  };

  gStyle->SetOptStat(0);
  TCanvas * matCanv = new TCanvas("matCanv","matCanv",1000,1000);
  orthMatrix->Draw("colztext");

};


Double_t Azimod(Int_t twist_, Int_t m_, Float_t phiH_, Float_t phiR_) {
  Double_t retval_;
  Double_t norm_ = 1 / ( TMath::Sqrt(2) * PI );
  if(twist_==0) return 1 / ( 2*PI ); // (constant)
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

  Bool_t disableLegendre = 1;

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
  return retval_ * norm_;
};
