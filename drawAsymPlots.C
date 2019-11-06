// draw asymmetry plots, with colors and formatting

void drawAsymPlots(
  TString infileN = "forMC/rec.02.root"
) {

  TFile * infile = new TFile(infileN,"READ");

  enum setting_enum {kOne,kMulti,kExtra};
  int setting = kMulti;

  TString pngName = infileN;
  TString pdfName = infileN;
  pngName(TRegexp("root$"))="";
  pdfName(TRegexp("root$"))="";

  TIter nextKey(gDirectory->GetListOfKeys());
  TString keyname;
  TCanvas * canv;
  TCanvas * canv2;
  TMultiGraph * mg;
  TLatex * text1 = new TLatex(0.14,0.85,"CLAS Preliminary");
  TLatex * text2 = new TLatex(0.14,0.8,"e^{-}p #rightarrow e^{-}#pi^{+}#pi^{-}X");

  gStyle->SetLegendTextSize(0.05);
  gStyle->SetLabelFont(22,"xyz");
  gStyle->SetTitleFont(22,"xyz");
  gStyle->SetTitleFont(22,"");

  //TLegend * leg = new TLegend(0.55, 0.1, 0.9, 0.32);
  TLegend * leg = new TLegend(0.6, 0.71, 0.93, 0.92);

  // hard coded plot nameing / coloring
  const Int_t numAmps = 4;
  TString legName[numAmps] = {
    "A_{R}sin #phi_{R}            ",
    "A_{hR}sin(#phi_{h}-#phi_{R})",
    "A_{h}sin #phi_{h}",
    "A_{hR}sin(#phi_{h}+#phi_{R})"
  };
  unsigned int color[numAmps] = {
    kGreen+1,
    kRed,
    kBlue,
    kMagenta
  };
  unsigned int markerStyle[numAmps] = {
    kFullTriangleUp,
    kFullCircle,
    kFullTriangleDown,
    kFullCross
  };
  unsigned int lineStyle[numAmps] = { 2, 1, 3, 4 };
  Bool_t legFilled = false;
  Int_t l;
  text1->SetNDC(1);
  text2->SetNDC(1);
  TString grT;
  Float_t xmin,xmax;
  TLine * zline;

  while(TKey * k = (TKey*) nextKey()) {
    keyname = TString(k->GetName());

    //if( keyname.Contains(TRegexp("^RF_.*_kindep_")) ) {
    if( keyname.Contains(TRegexp("^multiGr")) ) {
      mg = (TMultiGraph*) k->ReadObj();

      TIter nextGraph(mg->GetListOfGraphs());
      l=0;
      while(TGraph * g = (TGraph*) nextGraph()) {
        g->SetMarkerSize(0.6);
        if(setting==kMulti) {
          g->SetMarkerColor(color[l]);
          g->SetLineColor(color[l]);
          g->SetMarkerColor(color[l]);
          g->SetLineStyle(lineStyle[l]);
          g->SetMarkerStyle(markerStyle[l]);
          //g->SetLineStyle(1);
          if(!legFilled) leg->AddEntry(g,l<numAmps?legName[l]:"","PLE");
        } else if(setting==kOne) {
          g->SetMarkerStyle(kFullCircle);
          g->SetMarkerColor(kBlack);
          g->SetLineStyle(1);
        } else if(setting==kExtra) {
          if(l==3) { 
            g->SetMarkerStyle(kCircle); 
            g->SetLineColor(kBlack); 
          };
        };
        l++;
      };
      legFilled = true; // stop multidim binning from filling legend multiple times
      
      grT = TString(mg->GetTitle());
      //grT(TRegexp("_{m.l.m.}"))="";
      grT(TRegexp("A_{LU}.*vs")) = "A_{LU} vs";
      grT(TRegexp(" ::.*77)")) = ", for M_{h}<M_{#rho}";
      grT(TRegexp(" ::.*3.00)")) = ", for M_{h}>M_{#rho}";
      mg->SetTitle(grT);


      canv = new TCanvas("canv","canv",300,300);
      //canv->SetGrid(0,1);
      
      if(setting==kOne) mg->Draw("APE");
      else if(setting==kMulti || setting==kExtra) mg->Draw("LAPE");

      xmin = mg->GetXaxis()->GetXmin();
      xmax = mg->GetXaxis()->GetXmax();
      zline = new TLine(xmin,0,xmax,0);
      zline->SetLineColor(kGray+2);
      zline->SetLineWidth(1);

      mg->GetYaxis()->SetRangeUser(-0.2,0.3);
      mg->GetXaxis()->SetLabelSize(0.045);
      mg->GetYaxis()->SetLabelSize(0.045);
      mg->GetXaxis()->SetTitleSize(0.05);
      mg->GetXaxis()->SetTitleOffset(0.85);
      zline->Draw();

      //mg->GetYaxis()->UnZoom(); mg->Draw("LAPE");

      if(setting!=kOne) leg->Draw();

      //text1->Draw();
      //text2->Draw();

      //canv->Print(TString(pngName+keyname+".png"),"png"); // DO NOT USE!
      canv->Print(TString(pdfName+keyname+".pdf"),"pdf"); // instead use imagemagick
                                                          // to convert pdf to png
    };

    /*
    if( keyname.Contains(TRegexp("^chindfCanv")) ||
        keyname.Contains(TRegexp("^rellumCanv"))
      ) {
      canv2 = (TCanvas*) k->ReadObj();
      canv2->Print(TString(pdfName+keyname+".pdf"),"pdf");
      //canv2->Print(TString(pngName+keyname+".png"),"png");
    };
    */

  };
};


