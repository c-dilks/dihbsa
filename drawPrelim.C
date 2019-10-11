void drawPrelim(
  TString infileN = "forPrelim/m2.0x34.i2.three.root"
) {

  TFile * infile = new TFile(infileN,"READ");

  TString pngName = infileN;
  TString pdfName = infileN;
  pngName(TRegexp("root$"))="";
  pdfName(TRegexp("root$"))="";

  TIter nextKey(gDirectory->GetListOfKeys());
  TString keyname;
  TCanvas * canv;
  TCanvas * canv2;
  TMultiGraph * gr;
  TLatex * text1 = new TLatex(0.14,0.85,"CLAS Preliminary");
  TLatex * text2 = new TLatex(0.14,0.8,"e^{-}p #rightarrow e^{-}#pi^{+}#pi^{-}X");
  text1->SetNDC(1);
  text2->SetNDC(1);
  TString grT;
  Float_t xmin,xmax;
  TLine * zline;

  while(TKey * k = (TKey*) nextKey()) {
    keyname = TString(k->GetName());

    //if( keyname.Contains(TRegexp("^RF_.*_kindep_")) ) {
    if( keyname.Contains(TRegexp("^multiGr")) ) {
      gr = (TMultiGraph*) k->ReadObj();

      TIter nextGraph(gr->GetListOfGraphs());
      while(TGraph * g = (TGraph*) nextGraph()) {
        g->SetMarkerSize(0.6);
      };
      
      grT = TString(gr->GetTitle());
      //grT(TRegexp("_{m.l.m.}"))="";
      grT(TRegexp("A_{LU}.*vs")) = "A_{LU} vs";
      gr->SetTitle(grT);

      canv = new TCanvas("canv","canv",300,300);
      //canv->SetGrid(0,1);
      gr->Draw("LAPE");

      xmin = gr->GetXaxis()->GetXmin();
      xmax = gr->GetXaxis()->GetXmax();
      zline = new TLine(xmin,0,xmax,0);
      zline->SetLineColor(kGray+2);
      zline->SetLineWidth(1);

      gr->GetYaxis()->SetRangeUser(-0.08,0.08);
      //gr->GetYaxis()->SetRangeUser(-0.3,0.3);
      //gr->GetYaxis()->SetRangeUser(-0.16,0.16);
      zline->Draw("same");
      //gr->GetYaxis()->UnZoom(); gr->Draw("LAPE");
      text1->Draw();
      text2->Draw();
      //canv->Print(TString(pngName+keyname+".png"),"png"); // DO NOT USE!
      canv->Print(TString(pdfName+keyname+".pdf"),"pdf"); // instead use imagemagick
                                                          // to convert pdf to png
    };

    if( keyname.Contains(TRegexp("^chindfCanv")) ||
        keyname.Contains(TRegexp("^rellumCanv"))
      ) {
      canv2 = (TCanvas*) k->ReadObj();
      canv2->Print(TString(pdfName+keyname+".pdf"),"pdf");
      canv2->Print(TString(pngName+keyname+".png"),"png");
    };

  };
};


