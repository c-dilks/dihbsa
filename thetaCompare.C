// compare theta definitions

R__LOAD_LIBRARY(DihBsa)

#include "Constants.h"
#include "EventTree.h"
#include "Tools.h"
#include "TLorentzVector.h"
#include "TString.h"
#include "TMath.h"
#include "Dihadron.h"

void thetaCompare(TString dir="outroot.fall18.one") {

  TString files = dir + "/*.root";
  /////////////////////////////////
  Int_t whichPair = EncodePairType(kPip,kPim);
  /////////////////////////////////

  EventTree * ev = new EventTree(files,whichPair);


  TFile * outfile = new TFile("theta.root","RECREATE");
  TTree * tr = new TTree("tr","tr");
  Dihadron * dih;
  Float_t theta,thetaAlt,thetaLI,zeta,Mh;
  tr->Branch("evnum",&(ev->evnum));
  tr->Branch("theta",&theta);
  tr->Branch("thetaAlt",&thetaAlt);
  tr->Branch("thetaLI",&thetaLI);
  tr->Branch("zeta",&zeta);
  tr->Branch("Mh",&Mh);
  
  for(int i=0; i<ev->ENT; i++) {
    ev->GetEvent(i);

    if(ev->Valid()) {

      dih = ev->GetDihadronObj();

      theta = dih->theta;
      thetaAlt = dih->thetaAlt;
      thetaLI = dih->thetaLI;
      zeta = dih->zeta;
      Mh = dih->Mh;
      tr->Fill();

      /*
      printf("%d: theta=%f  thetaAlt=%f  thetaLI=%f\n",
        ev->evnum, theta, thetaAlt, thetaLI);
      */

    };
  };

  tr->Write();
  outfile->Close();


};
