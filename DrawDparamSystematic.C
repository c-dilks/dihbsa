// draws plots for determining systematic from D_1 pp-wave
void DrawDparamSystematic(TString dirN = "spinroot") {

  // get list of "Dasym" files
  TObjArray * asymFileArr = new TObjArray();
  TSystemDirectory dir(dirN,dirN);
  TList * asymFiles = dir.GetListOfFiles();
  TString asymFileN;
  if(asymFiles) {
    TIter nextSF(asymFiles);
    TSystemFile * asymFile;
    while((asymFile=(TSystemFile*)nextSF())) {
      asymFileN = asymFile->GetName();
      if(!asymFile->IsDirectory() && 
          asymFileN.Contains("Dasym") &&
          asymFileN.EndsWith(".root")) {
        asymFileArr->AddLast(new TFile(TString(dirN+"/"+asymFileN),"READ"));
      };
    };
  };

  
  // arrays of graphs; one entry for one asymmetry amplitude
  TObjArray * asymGrArr = new TObjArray(); // for "control" asymmetry graphs (Dparam=0)
  const Int_t NBINS = 10;
  TObjArray * DGrArr[NBINS]; // for asym vs. Dparam, where array index is for each kin. bin
  int b;
  for(b=0; b<NBINS; b++) DGrArr[b] = new TObjArray();

  TString DparamStr;
  Float_t Dparam;


  // loop through Dasym files
  TFile * infile;
  TIter next(asymFileArr);
  TKey * key;
  TString keyname;
  Bool_t controlFile;
  Bool_t first = true;
  TGraphErrors * gr;
  TGraphErrors * dgr;
  Int_t amp;
  TString ampStr;
  while((infile=(TFile*)next())) {

    // focus the file
    printf("read file %s\n",infile->GetName());
    infile->cd();

    // get Dparam
    DparamStr = infile->GetName();
    DparamStr(TRegexp("^.*_")) = "";
    DparamStr(TRegexp(".root$")) = "";
    Dparam = DparamStr.Atof();
    printf("-> Dparam = %f\n",Dparam);
    controlFile = fabs(Dparam)<0.001;

    // get asymmetry graphs
    TIter nextKey(gDirectory->GetListOfKeys());
    while((key=(TKey*)nextKey())) {
      keyname = TString(key->GetName());
      if(keyname.Contains(TRegexp("^RF_")) && 
         keyname.Contains("_kindep_") &&
         TString(key->GetClassName())=="TGraphErrors") {

        //printf(" %s %s\n",keyname.Data(),key->GetClassName());
        gr = (TGraphErrors*) key->ReadObj();
        if(controlFile) asymGrArr->AddLast(gr);

        // get the amplitude number
        TString ampStr = keyname;
        ampStr(TRegexp("^RF_A")) = "";
        ampStr(TRegexp("_.*$")) = "";
        amp = ampStr.Atoi();

        // instantiate asymVsDparam graphs
        if(first) {
          for(int i=0; i<gr->GetN(); i++) {
            dgr = new TGraphErrors();
            dgr->SetName(TString("asymVsDparam_"+TString::Itoa(i,10)));
            DGrArr[i]->AddLast(dgr);
          };
        };

        // fill asymVsDparam graphs
        for(int i=0; i<gr->GetN(); i++) {
          DGrArr[i]->At(amp)->SetPoint(...);
        };




      };
    };

    first = false;

  };

};


