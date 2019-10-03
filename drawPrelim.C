void drawPrelim(
  TString infileN = "forPrelim/m1.0x34.i1.root"
) {

  TFile * infile = new TFile(infileN,"READ");

  TString pngName = infileN;
  TString pdfName = infileN;
  pngName(TRegexp("root$"))="png";
  pdfName(TRegexp("root$"))="pdf";

  TIter nextKey(gDirectory->GetListOfKeys());
  TString keyname;
  TCanvas * canv;
  TGraphErrors * gr;
  TString grT;

  while(TKey * k = (TKey*) nextKey()) {
    keyname = TString(k->GetName());

    if( keyname.Contains(TRegexp("^RF_.*_kindep_")) ) {
      cout << "found: " << keyname << endl;
      gr = (TGraphErrors*) k->ReadObj();
      
      grT = TString(gr->GetTitle());
      grT(TRegexp("_{m.l.m.}"))="";
      gr->SetTitle(grT);

      canv = new TCanvas("canv","canv",800,800);
      canv->SetGrid(1,1);
      gr->Draw("APE");
      //gr->GetYaxis()->UnZoom();
      gr->Draw("APE");
      canv->Print(pngName,"png");
      //canv->Print(pdfName,"pdf");
    };

  };
};


