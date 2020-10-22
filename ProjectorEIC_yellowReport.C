// projections for yellow report
//
// - produces `dihadronPWprojection.png`
//   - compares two different minimum pion pT values
//   - to disable one or the other, prevent it from being
//     added to `mgr`
// - the vector file `dihadronPWprojection.svg` links to
//   the png `dihadronPWprojection.png`, and contains latex
//   overlays for the labels
// - latex overlays made with tex2img, font size 42
// - export svg file as png with dpi=100

void ProjectorEIC_yellowReport(
  TString infile0N="spinroot_5x41_300/asym_test1.root",
  TString infile1N="spinroot_5x41_100/asym_test1.root"
) {

  // read asymmetry graphs
  int f;
  TFile * infile[2];
  infile[0] = new TFile(infile0N,"READ");
  infile[1] = new TFile(infile1N,"READ");
  TListIter nextKey(infile[0]->GetListOfKeys());
  TString keyname;
  TObjArray * asymArr[2];
  for(f=0; f<2; f++) asymArr[f] = new TObjArray();
  while(TKey * key = (TKey*) nextKey()) {
    keyname = TString(key->GetName());
    // read asymmetry graph
    if(keyname.Contains(TRegexp("^kindepMA")) &&
           !keyname.Contains("Canv")) {
      for(f=0; f<2; f++) {
        asymArr[f]->AddLast((TGraphAsymmErrors*)infile[f]->Get(keyname));
      };
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
  printf("number of generated events: %f\n",evCountGen);
  printf("cross section: %f nb\n",crossSection*1e9);
  printf("generated luminosity: %f fb^-1\n",luminosityGen/1e15);
  printf("scale factor to %f fb^-1: %f\n",luminosityProj/1e15,scaleFactor);
  //printf("number of events for %f fb^-1: %f\n",luminosity,evCountProj);
  //printf("yield scale factor: %f\n",evCountProj/evCountGen);
  printf("scale uncertainty by: %f\n",errScale);


  // project asymmetry graphs
  TMultiGraph * mgr;
  TGraphAsymmErrors * gr;
  TString grT;
  Double_t x,y,ex,ey;
  TLine * zero;
  TCanvas * canv = new TCanvas("canv","canv",3600,1800);
  Int_t pad = 1;
  Int_t amp;
  char noop[32];
  TString mod,title,eigen,ffpol,super;
  //TLatex * eigenTex;
  //TLatex * ffpolTex;
  canv->Divide(3,3);
  for(int e=0; e<asymArr[0]->GetEntries(); e++) {
    for(f=0; f<2; f++) {
      gr = (TGraphAsymmErrors*) asymArr[f]->At(e);
      gr->SetMarkerColor(kBlack);
      gr->SetLineColor(f==0?kBlue:kRed);
      gr->SetFillColor(f==0?kBlue:kRed);
      gr->SetMarkerSize(0);
      gr->SetMarkerStyle(kFullCircle);
      //gr->SetLineWidth(f==0?2:5); // pdf
      gr->SetLineWidth(f==0?10:20); // png

      sscanf(gr->GetName(),"kindepMA_A%d_%s",&amp,noop);
      cout << gr->GetName() << " " << amp << endl;
      switch(amp) {
        case 0:
          mod="sin(#phi_{h}+#phi_{S})";
          eigen = "0,0";
          ffpol = "OO";
          super = "#perpss+pp";
          break;
        case 1:
          mod="cos(#theta) sin(#phi_{h}+#phi_{S})";
          eigen = "1,0";
          ffpol = "OL";
          super = "#perp";
          break;
        case 2:
          mod="sin(#theta) sin(#phi_{R}+#phi_{S})";
          eigen = "1,1";
          ffpol = "OT";
          super = "#angle";
          break;
        case 3:
          mod="sin(#theta) sin(2#phi_{h}-#phi_{R}+#phi_{S})";
          eigen = "1,-1";
          ffpol = "OT";
          super = "#perp";
          break;
        case 4:
          mod="1/2 (3cos^{2}#theta-1) sin(#phi_{h}+#phi_{S})";
          eigen = "2,0";
          ffpol = "LL";
          super = "#perp";
          break;
        case 5:
          mod="sin(2#theta) sin(#phi_{R}+#phi_{S})";
          eigen = "2,1";
          ffpol = "LT";
          super = "#angle";
          break;
        case 6:
          mod="sin(2#theta) sin(2#phi_{h}-#phi_{R}+#phi_{S})";
          eigen = "2,-1";
          ffpol = "LT";
          super = "#perp";
          break;
        case 7:
          mod="sin^{2}(#theta) sin(-#phi_{h}+2#phi_{R}+#phi_{S})";
          eigen = "2,2";
          ffpol = "TT";
          super = "#angle";
          break;
        case 8:
          mod="sin^{2}(#theta) sin(3#phi_{h}-2#phi_{R}+#phi_{S})";
          eigen = "2,-2";
          ffpol = "TT";
          super = "#perp";
          break;
      };
      title = gr->GetTitle();
      mod = "A_{UT}^{"+mod+"} vs.";
      eigen = "|"+eigen+">";
      ffpol = "h_{1}#otimesH_{"+ffpol+"}^{"+super+"}";
      title(TRegexp("A_{UT}.* vs.")) = mod;
      //eigenTex = new TLatex(0.75,0.25,eigen);
      //ffpolTex = new TLatex(0.75,0.15,ffpol);
      //eigenTex->SetNDC(1);
      //ffpolTex->SetNDC(1);
      //eigenTex->SetTextSize(0.08);
      //ffpolTex->SetTextSize(0.08);

      gStyle->SetTitleBorderSize(4);
      gStyle->SetTitleSize(0.08,"main");
      gStyle->SetTitleStyle(1001);
      gStyle->SetTitleFillColor(kWhite);

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

      if(f==0) {
        mgr = new TMultiGraph();
        mgr->SetTitle(title);
        mgr->GetYaxis()->SetTitle("A_{UT}");
        mgr->GetXaxis()->SetTitle("x");
        //mgr->GetYaxis()->SetNoExponent(1);
        mgr->GetXaxis()->SetTitleOffset(0.85);
        mgr->GetYaxis()->SetTitleOffset(0.85);
        mgr->GetXaxis()->SetTitleSize(0.06);
        mgr->GetYaxis()->SetTitleSize(0.06);
        mgr->GetXaxis()->SetLabelSize(0.06);
        mgr->GetYaxis()->SetLabelSize(0.06);
      };
      gr->GetXaxis()->SetLabelSize(0.06);
      gr->GetYaxis()->SetLabelSize(0.06);

      if(f==1) mgr->Add(gr); // 100 MeV only
      //mgr->Add(gr); // 100 MeV only

    };

    //zero = new TLine(gr->GetXaxis()->GetXmin(),0,gr->GetXaxis()->GetXmax(),0);
    zero = new TLine(0,0,0.14,0);
    //zero->SetLineStyle(2);
    zero->SetLineWidth(1);

    canv->cd(pad);
    canv->GetPad(pad)->SetGrid(1,1);
    //mgr->Draw("A3"); // error band
    mgr->Draw("APEZ"); // error bars
    mgr->GetXaxis()->SetLimits(0,0.14); // for 5x41
    mgr->GetYaxis()->SetRangeUser(-9e-4,9e-4); // for 5x41 test 1
    //mgr->GetYaxis()->SetRangeUser(-2e-4,2e-4); // for 5x41 tests 2 & 3
    //mgr->GetYaxis()->SetRangeUser(-5e-3,5e-3); // for 18x275
    zero->Draw();

    //eigenTex->Draw();
    //ffpolTex->Draw();
    pad++;
  };

  canv->Print("dihadronPWprojection.png","png");
};
    
