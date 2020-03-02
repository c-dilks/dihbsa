// takes Timothy's grids, from mass_*.txt, converts them into histograms, and inserts
// them into a cat.root file, overwriting the ones that are already in there
// - this is used to read Timothy's grid into a cat.root file, that I can easily feed to
//   asymFit.exe
// - PROCEDURE:
//   - make sure we have a `cat.root` file in the `spinroot` directory
//   - cp my `spinroot` directory here as `spinroot.chrisGrids`
//   - cp the `cat.root` file to `spinroot.timGrids` as well; the `cat.root` file here
//     will have its grids updated, according to the text-file grids from Timothy, stored
//     in the `grids` subdirectory

Int_t nx,ny;
Double_t xl,xu,yl,yu;

TH2D * mkHisto(TString fN, TString suf) {
  TTree * t = new TTree();
  t->ReadFile(fN,"phiR/D:phiH/D:c/D");
  Double_t phiH,phiR,c;
  t->SetBranchAddress("phiH",&phiH);
  t->SetBranchAddress("phiR",&phiR);
  t->SetBranchAddress("c",&c);
  Int_t nbins = 8;
  TString hN = "stream_cat_aziDist_"+suf;
  TString hT = "grid: "+fN+";#phi_{R};#phi_{h}";
  //TH2D * h = new TH2D(hN,hT,nbins,0,2*TMath::Pi(),nbins,0,2*TMath::Pi());
  TH2D * h = new TH2D(hN,hT,nx,xl,xu,ny,yl,yu);
  for(int i=0; i<t->GetEntries(); i++) {
    t->GetEntry(i);
    //h->Fill(phiR,phiH,c);
    for(int k=0; k<c; k++) h->Fill(phiR,phiH);
  };
  return h;
};

//////////////////////
void embedGrid() {

  TFile * oldfile = new TFile("spinroot.chrisGrids/cat.root","READ");
  TH2D * oldhist = (TH2D*) oldfile->Get("A_M0/stream_cat_aziDist_P_M0");
  nx = oldhist->GetXaxis()->GetNbins();
  xu = oldhist->GetXaxis()->GetXmax();
  xl = oldhist->GetXaxis()->GetXmin();
  ny = oldhist->GetYaxis()->GetNbins();
  yu = oldhist->GetYaxis()->GetXmax();
  yl = oldhist->GetYaxis()->GetXmin();
  oldfile->Close();


  const Int_t NH = 2;
  const Int_t NM = 2;
  TH2D * hist[NM][NH]; // [mass bin] [hel+-]
  ///*
  hist[0][0] = mkHisto("grids/mass_posHel_lowerBin.txt","P_M0");
  hist[0][1] = mkHisto("grids/mass_negHel_lowerBin.txt","M_M0");
  hist[1][0] = mkHisto("grids/mass_posHel_upperBin.txt","P_M1");
  hist[1][1] = mkHisto("grids/mass_negHel_upperBin.txt","M_M1");
  //*/
  /*
  hist[0][0] = mkHisto("printGrid/timothy_m0_hp.txt","P_M0");
  hist[0][1] = mkHisto("printGrid/timothy_m0_hm.txt","M_M0");
  hist[1][0] = mkHisto("printGrid/timothy_m1_hp.txt","P_M1");
  hist[1][1] = mkHisto("printGrid/timothy_m1_hm.txt","M_M1");
  */
  
  TFile * infile = new TFile("spinroot.timGrids/cat.root","UPDATE");
  TString writeN,dirN;
  for(int m=0; m<NM; m++) {
    for(int h=0; h<NH; h++) {
      dirN = Form("/A_M%d",m);
      infile->cd(dirN);
      writeN = hist[m][h]->GetName();
      hist[m][h]->Write(writeN,TObject::kOverwrite);
    };
  };
  infile->Close();
};
