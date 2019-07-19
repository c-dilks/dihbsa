void PrintEnergies() {
  TFile * infile = new TFile("../outroot/out_clasdispr.00.e10.600.emn0.75tmn.09.xs80.53nb.dis.0000.nrad.dat.evio.root","READ");
  TTree * t = (TTree*) infile->Get("tree");
  Int_t evnum,pairType;
  Float_t hadE[2];
  t->SetBranchAddress("evnum",&evnum);
  t->SetBranchAddress("hadE",hadE);
  t->SetBranchAddress("pairType",&pairType);

  for(int x=0; x<t->GetEntries(); x++) {
    t->GetEntry(x);
    if(pairType==0x34) {   
      printf("%d %f %f\n",evnum,hadE[0],hadE[1]);
    };
  };
};

