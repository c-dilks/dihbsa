void ProjectorEIC(TString infileN="spinroot_5x41_100/asym_test1.root") {

  // read asymmetry graphs
  TFile * infile = new TFile(infileN,"READ");
  TListIter nextKey(gDirectory->GetListOfKeys());
  TString keyname;
  Double_t evCountRec = 0;
  TObjArray * asymArr = new TObjArray();
  while(TKey * key = (TKey*) nextKey()) {
    keyname = TString(key->GetName());
    // get number of events (actually, dihadrons) reconstructed from the MC
    if(keyname.Contains(TRegexp("^ivFullDist"))) {
      evCountRec += ((TH1D*)key->ReadObj())->GetEntries();
    }
    // read asymmetry graph
    else if(keyname.Contains(TRegexp("^kindepMA")) &&
           !keyname.Contains("Canv")) {
      asymArr->AddLast((TGraphAsymmErrors*)key->ReadObj());
    };
  };


  // error scale factor
  Double_t evCountGen = 1e6;
  Double_t crossSection = 0.36241800884389658 * 1e-6; // [barns] // for 5x41
  Double_t luminosityGen = evCountGen / crossSection; // [barns^-1]
  Double_t luminosityProj = 10.0/1e-15; // [barns^-1]
  Double_t scaleFactor = luminosityProj / luminosityGen;
  //Double_t evCountProj = crossSection * 1e-6 * luminosity / 1e-15;
  Double_t errScale = 1.0 / TMath::Sqrt(scaleFactor);
  printf("number of reconstructed dihadrons: %f\n",evCountRec);
  printf("number of generated events: %f\n",evCountGen);
  printf("cross section: %f nb\n",crossSection*1e9);
  printf("generated luminosity: %f fb^-1\n",luminosityGen/1e15);
  printf("scale factor to %f fb^-1: %f\n",luminosityProj/1e15,scaleFactor);
  //printf("number of events for %f fb^-1: %f\n",luminosity,evCountProj);
  //printf("yield scale factor: %f\n",evCountProj/evCountGen);
  printf("scale uncertainty by: %f\n",errScale);


  // project asymmetry graphs
  TObjArrayIter nextGr(asymArr);
  TString grT;
  Double_t x,y,ex,ey;
  TLine * zero;
  Int_t ngraphs = asymArr->GetEntries();
  TCanvas * canv = new TCanvas("canv","canv",3600,1800);
  Int_t pad = 1;
  Int_t amp;
  char noop[32];
  TString mod,title;
  canv->Divide(3,3);
  while(TGraphAsymmErrors * gr = (TGraphAsymmErrors*) nextGr()) {
    gr->SetMarkerColor(kBlack);
    gr->SetLineColor(kBlue);
    gr->SetFillColor(kBlue);
    gr->SetMarkerSize(0);
    gr->SetMarkerStyle(kFullCircle);
    gr->SetLineWidth(4);

    sscanf(gr->GetName(),"kindepMA_A%d_%s",&amp,noop);
    cout << gr->GetName() << " " << amp << endl;
    switch(amp) {
      case 0: mod="sin(#phi_{h}+#phi_{S})"; break;
      case 1: mod="cos(#theta) sin(#phi_{h}+#phi_{S})"; break;
      case 2: mod="sin(#theta) sin(#phi_{R}+#phi_{S})"; break;
      case 3: mod="sin(#theta) sin(2#phi_{h}-#phi_{R}+#phi_{S})"; break;
      case 4: mod="1/2 (3cos^{2}#theta-1) sin(#phi_{h}+#phi_{S})"; break;
      case 5: mod="sin(2#theta) sin(#phi_{R}+#phi_{S})"; break;
      case 6: mod="sin(2#theta) sin(2#phi_{h}-#phi_{R}+#phi_{S})"; break;
      case 7: mod="sin^{2}(#theta) sin(-#phi_{h}+2#phi_{R}+#phi_{S})"; break;
      case 8: mod="sin^{2}(#theta) sin(3#phi_{h}-2#phi_{R}+#phi_{S})"; break;
    };
    title = gr->GetTitle();
    mod = "A_{UT}^{"+mod+"} vs.";
    title(TRegexp("A_{UT}.* vs.")) = mod;
    gr->SetTitle(title);
    gr->GetYaxis()->SetTitle("A_{UT}");

    gStyle->SetTitleBorderSize(4);
    gStyle->SetTitleSize(0.08,"main");
    gStyle->SetTitleStyle(1001);
    gStyle->SetTitleFillColor(kWhite);
    //gr->GetYaxis()->SetNoExponent(1);
    gr->GetXaxis()->SetTitleOffset(0.85);
    gr->GetYaxis()->SetTitleOffset(0.85);
    gr->GetXaxis()->SetTitleSize(0.06);
    gr->GetYaxis()->SetTitleSize(0.06);
    gr->GetXaxis()->SetLabelSize(0.06);
    gr->GetYaxis()->SetLabelSize(0.06);

    for(int i=0; i<gr->GetN(); i++) {
      gr->GetPoint(i,x,y);
      ex = gr->GetErrorX(i); // parabolic error from HESSE
      ey = gr->GetErrorY(i);

      // divide polarization
      ey /= 0.7;

      // divide depolarization ratio <B/A>
      switch(i) {
	      case 0: ey /= 0.60614853; break;
	      case 1: ey /= 0.84335763; break;
	      case 2: ey /= 0.92392810; break;
	      case 3: ey /= 0.96240604; break;
	      case 4: ey /= 0.97970662; break;
        default: fprintf(stderr,"ERROR: unknown depol factor\n");
      };


      gr->SetPoint(i,x,0);
      gr->SetPointError(i,ex,ex,ey*errScale,ey*errScale); // (TGraphAsymmErrors)
    };

    //zero = new TLine(gr->GetXaxis()->GetXmin(),0,gr->GetXaxis()->GetXmax(),0);
    zero = new TLine(0,0,0.14,0);
    //zero->SetLineStyle(2);
    zero->SetLineWidth(1);

    canv->cd(pad);
    canv->GetPad(pad)->SetGrid(1,1);
    //gr->Draw("A3"); // error band
    gr->Draw("APEZ"); // error bars
    gr->GetXaxis()->SetLimits(0,0.14); // for 5x41
    gr->GetYaxis()->SetRangeUser(-7e-4,7e-4); // for 5x41 test 1
    //gr->GetYaxis()->SetRangeUser(-2e-4,2e-4); // for 5x41 tests 2 & 3
    //gr->GetYaxis()->SetRangeUser(-5e-3,5e-3); // for 18x275
    zero->Draw();
    pad++;
  };

  canv->Print("dihadronPWprojection.png","png");
};
    
