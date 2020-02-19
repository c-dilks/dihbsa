// draws grids from two different cat.root files, and fits them
// - a 'grid' is the yield binned in (phiH,phiR)
// - to make a grid, you need to use buildSpinroot.exe with the -b option

void drawGrid(TString catFileC="spinroot.chrisGrids/cat.root",
              TString catFileT="spinroot.timGrids/cat.root") {
  enum e {chris,tim};
  TFile * infile[2];
  infile[chris] = new TFile(catFileC,"READ");
  infile[tim] = new TFile(catFileT,"READ");

  const Int_t NH = 2;
  const Int_t NM = 2;
  TH2D * hist[2][NM][NH]; // [chris,tim] [mass bin] [hel-+]
  TString histN,dirN;

  for(int n=0; n<2; n++) {
    infile[n]->cd();
    for(int m=0; m<NM; m++) {
      for(int h=0; h<NH; h++) {
        dirN = Form("/A_M%d",m);
        histN = Form("/A_M%d/stream_cat_aziDist_%s_M%d",m,h==0?"M":"P",m);
        hist[n][m][h] = (TH2D*) infile[n]->Get(histN);
        printf("%p\n",(void*)hist[n][m][h]);
      };
    };
  };

  TCanvas * canv[2];
  canv[chris] = new TCanvas("canvChris","canvChris",1000,1000);
  canv[tim] = new TCanvas("canvTimothy","canvTimothy",1000,1000);
  gStyle->SetOptStat(0);
  
  for(int n=0; n<2; n++) {
    canv[n]->Divide(2,2);
    for(int m=0; m<NM; m++) {
      for(int h=0; h<NH; h++) {
        canv[n]->cd(2*h+m+1);
        hist[n][m][h]->Draw("colztext");
      };
    };
  };

  TH2D * comp[NM][NH]; // [mass bin] [hel-+]
  TCanvas * canvComp = new TCanvas("canvComp","canvComp",1000,1000);
  canvComp->Divide(2,2);
  TString compN,compT;
  for(int m=0; m<NM; m++) {
    for(int h=0; h<NH; h++) {
      compN = Form("comp_%d_%d",m,h);
      comp[m][h] = (TH2D*) hist[chris][m][h]->Clone(compN);
      //comp[m][h]->SetName(compN);
      compT = "chris/timothy :: " + TString(hist[chris][m][h]->GetTitle());
      comp[m][h]->SetTitle(compT);
      comp[m][h]->Divide(hist[tim][m][h]); // ignore complates about different axis lims
      canvComp->cd(2*h+m+1);
      comp[m][h]->SetMinimum(1);
      comp[m][h]->SetMaximum(1.2);
      comp[m][h]->Draw("colztext");
    };
  };

  canv[chris]->Print("gridChris.png","png");
  canv[tim]->Print("gridTimothy.png","png");
  canvComp->Print("gridRatio.png","png");

};
