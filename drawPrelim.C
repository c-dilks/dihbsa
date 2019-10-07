void drawPrelim(
  TString infileN = "forPrelim/m1.0x34.i1.root"
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
  TString grT;

  while(TKey * k = (TKey*) nextKey()) {
    keyname = TString(k->GetName());

    //if( keyname.Contains(TRegexp("^RF_.*_kindep_")) ) {
    if( keyname.Contains(TRegexp("^multiGr")) ) {
      gr = (TMultiGraph*) k->ReadObj();
      
      grT = TString(gr->GetTitle());
      //grT(TRegexp("_{m.l.m.}"))="";
      grT(TRegexp("A_{LU}.*vs")) = "A_{LU} vs";
      gr->SetTitle(grT);

      canv = new TCanvas("canv","canv",800,800);
      canv->SetGrid(0,1);
      gr->Draw("LAPE");
      gr->GetYaxis()->SetRangeUser(-0.11,0.11);
      gr->Draw("LAPE");
      //gr->GetYaxis()->UnZoom(); gr->Draw("LAPE");
      canv->Print(TString(pngName+keyname+".png"),"png");
      canv->Print(TString(pdfName+keyname+".pdf"),"pdf");
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


