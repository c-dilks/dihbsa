void mkTree(TString f="missingEvents.dat") {
  TTree * t = new TTree("t","t");
  TString branches = "evnum/I:hadOrder/I";
  branches += ":pipP/F:pipPt/F:pipEta/F:pipPhi/F";
  branches += ":pimP/F:pimPt/F:pimEta/F:pimPhi/F";
  branches += ":eleP/F:elePt/F:eleEta/F:elePhi/F";
  branches += ":Q2/F:W/F:x/F:y/F:Mh/F:xF/F:PhPerp/F:theta/F:PhiH/F:PhiR/F";
  t->ReadFile(f,branches);
};
