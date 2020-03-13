// timothy has ~660 more events than I do; this looks in the outroot file and prints
// their data

void PrintMyMissingEvents() {
  TFile * myFile = new TFile("../outroot/skim5_5051.hipo.root","READ");
  TTree * myTree = (TTree*) myFile->Get("tree");
  TTree * missTree = new TTree();
  missTree->ReadFile("missingEvents.dat","evnum/I:pipP/F:pimP/F");
  Int_t evnum;
  missTree->SetBranchAddress("evnum",&evnum);
  TString cut;
  for(int i=0; i<missTree->GetEntries(); i++) {
    missTree->GetEvent(i);
    cut = Form("evnum==%d && pairType==0x34",evnum);
    myTree->Scan("evnum:hadP[0]:hadP[1]:Mmiss:Zpair",cut);
  };
};

