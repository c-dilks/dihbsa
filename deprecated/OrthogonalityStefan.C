// for studying Stefan's asymmetry denominator modulations' mutual orthogonality

Double_t Azimod(Int_t mode_, Float_t phi_);

void OrthogonalityStefan(TString infileN="ortho.root") {

  // open data hist
  TFile * infile = new TFile(infileN,"READ");
  TH1D * dataDist = (TH1D*) infile->Get("d1");
  dataDist->SetTitle("data distribution");

  int f,g,b;

  // get number of bins and bin widths
  Int_t nbinsPhi = dataDist->GetNbinsX();
  Float_t minPhi = dataDist->GetXaxis()->GetXmin();
  Float_t maxPhi = dataDist->GetXaxis()->GetXmax();


  // |<fg>| matrix
  const Int_t NF = 3;
  TH2D * orthMatrix = new TH2D("orthMatrix","<fg> matrix",NF,0,NF,NF,0,NF);
  TString funcT[NF] = {
    "const",
    "cos#phi",
    "cos2#phi"
  };
  for(f=0; f<NF; f++) { 
    orthMatrix->GetXaxis()->SetBinLabel(f+1,funcT[f]);
    orthMatrix->GetYaxis()->SetBinLabel(f+1,funcT[f]);
  };


  // calculate <fg>
  TH1D * intDist[NF][NF];
  Float_t phi;
  Double_t dataWeight,aziF,aziG,product,integral;
  TString intDistN,intDistT;
  Double_t weightedNorm[NF];
  Double_t valF,valG;
  for(f=0; f<NF; f++) {
    for(g=0; g<NF; g++) {

      // generate histos, which hold the product of data*f*g
      intDistN = Form("intDist_f%d_g%d",f,g);
      intDistT = "f: " + funcT[f] + ",  g: " + funcT[g];
      intDist[f][g] = new TH1D(intDistN,intDistT,nbinsPhi,minPhi,maxPhi);

      // loop over bins
      for(b=1; b<=nbinsPhi; b++) {
        phi = intDist[f][g]->GetXaxis()->GetBinCenter(b);
        aziF = Azimod(f,phi);
        aziG = Azimod(g,phi);
        dataWeight = dataDist->GetBinContent(b);
        //dataWeight = 1;
        product = dataWeight * aziF * aziG;
        intDist[f][g]->SetBinContent(b,product);
      };

      // computed weighted normalization for function f as 1/sqrt(<ff>)
      if(f==g) {
        integral = intDist[f][g]->Integral("width");
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
      //integral = TMath::Abs( intDist[f][g]->Integral("width") );
      integral = intDist[f][g]->Integral("width");
      orthMatrix->SetBinContent(f+1,g+1,integral);
    };
  };


  // draw
  gStyle->SetOptStat(0);
  //gStyle->SetPalette(kBird);
  gStyle->SetPaintTextFormat(".3f");
  TCanvas * matCanv = new TCanvas("matCanv","matCanv",1000,1000);
  orthMatrix->SetMinimum(-1);
  orthMatrix->SetMaximum(1);
  orthMatrix->Draw("boxtext");


  // print <fg> matrix (mathematica syntax)
  printf("A={\n");
  for(f=0; f<NF; f++) {
    printf("{");
    for(g=0; g<NF; g++) 
      printf("%.3f%s",orthMatrix->GetBinContent(f+1,g+1),g+1==NF?"":",");
    printf("}%s\n",f+1==NF?"":",");
  };
  printf("};\n\n");

  
  // draw intDist
  TCanvas * canvInt = new TCanvas("canvInt","canvInt",1000,1000);
  canvInt->Divide(NF,NF);
  Int_t pad;
  for(f=0; f<NF; f++) {
    for(g=0; g<NF; g++) {
      pad = NF*(NF-f-1)+g+1;
      canvInt->GetPad(pad)->SetGrid(0,1);
      canvInt->cd(pad);
      intDist[f][g]->Draw("hist");
    };
  };



};

Double_t Azimod(Int_t mode_, Float_t phi_) {
  if(mode_==0) return 1;
  else if(mode_==1) return TMath::Cos(phi_);
  else if(mode_==2) return TMath::Cos(2*phi_);
  else return -10000;
};
