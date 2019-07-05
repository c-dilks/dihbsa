// for printing kinematics to cross check with others

R__LOAD_LIBRARY(../src/DihBsa)

#include "Constants.h"
#include "Tools.h"
#include "EventTree.h"
#include <map>
#include <math.h>


// hash function
Float_t Digitize(Float_t val) { 
  Int_t digits = 1; // number of decimal digits to keep
  return roundf(val*pow(10,digits));
};
Int_t HashTim(Float_t Q2_, Float_t W_) { 
  return 1000 * Digitize(Q2_) + Digitize(W_);
};
Int_t HashHarut(Float_t pipE_, Float_t pimE_) {
  return 1000 * Digitize(pipE_) + Digitize(pimE_);
};


// printing methods
void PrintCompare(TString name, Float_t val, Float_t xval);


void CrossChecker(TString indir="../outroot") {

  ///////////////////////////////
  enum xenum { kTim, kHarut, kSimple };
  Int_t WHICH_XCHECK = kHarut;
  ///////////////////////////////

  Int_t whichPair = EncodePairType(kPip,kPim);
  EventTree * ev = new EventTree(TString(indir+"/*.root"),whichPair); 


  TFile * xfile;
  TTree * xtree = new TTree();
  TString xstr;

  Int_t evnum;
  Float_t Q2,W,x,y,Mh,pT,xF,theta,PhiH,PhiR;
  enum pippim_enum { kpip, kpim };
  Float_t hadE[2];
  Float_t hadPt[2];
  Float_t hadTheta[2];
  Float_t hadPhi[2];
  Float_t PhPt;

  Float_t valTheta;


  if(WHICH_XCHECK == kTim) {
    gROOT->ProcessLine(".! python formatTimFile.py");
    xstr = "evnum/I:Q2/F:W/F:x/F:y/F:Mh/F:pT/F:xF/F:theta/F:PhiR/F:PhiH/F";
    xtree->ReadFile("xtree.dat",xstr);
    xtree->SetBranchAddress("evnum",&evnum);
    xtree->SetBranchAddress("Q2",&Q2);
    xtree->SetBranchAddress("W",&W);
    xtree->SetBranchAddress("x",&x);
    xtree->SetBranchAddress("y",&y);
    xtree->SetBranchAddress("Mh",&Mh);
    xtree->SetBranchAddress("pT",&pT);
    xtree->SetBranchAddress("xF",&xF);
    xtree->SetBranchAddress("theta",&theta);
    xtree->SetBranchAddress("PhiH",&PhiH);
    xtree->SetBranchAddress("PhiR",&PhiR);
  }
  else if(WHICH_XCHECK == kHarut) {
    gROOT->ProcessLine(".! python formatHarutFile.py");
    xstr = "evnum/I";
    xstr += ":pipE/F:pipTheta/F:pipPhi/F:pipPt/F";
    xstr += ":pimE/F:pimTheta/F:pimPhi/F:pimPt/F";
    xstr += ":PhPt/F";
    xtree->ReadFile("xtree.dat",xstr);

    xtree->SetBranchAddress("evnum",&evnum);

    xtree->SetBranchAddress("pipE",&hadE[kpip]);
    xtree->SetBranchAddress("pipPt",&hadPt[kpip]);
    xtree->SetBranchAddress("pipTheta",&hadTheta[kpip]);
    xtree->SetBranchAddress("pipPhi",&hadPhi[kpip]);

    xtree->SetBranchAddress("pimE",&hadE[kpim]);
    xtree->SetBranchAddress("pimPt",&hadPt[kpim]);
    xtree->SetBranchAddress("pimTheta",&hadTheta[kpim]);
    xtree->SetBranchAddress("pimPhi",&hadPhi[kpim]);

    xtree->SetBranchAddress("PhPt",&PhPt);
  }
  else if(WHICH_XCHECK == kSimple) {
    xfile = new TFile("../simpleTree.root","READ");
    xtree = (TTree*) xfile->Get("tree");
    xtree->SetBranchAddress("evnum",&evnum);
    xtree->SetBranchAddress("hadE",hadE);
    xtree->SetBranchAddress("hadPt",hadPt);
  }
  else {
    fprintf(stderr,"ERROR: unknown WHICH_XCHECK\n");
    exit(0);
  };


  //xtree->Scan("*");


  // build hash table
  Float_t hashVal;
  std::map<Float_t,Int_t> hashMap;
  std::map<Float_t,Int_t>::iterator hashIter;
  for(int xi=0; xi<xtree->GetEntries(); xi++) {
    xtree->GetEntry(xi);

    switch(WHICH_XCHECK) {
      case kTim: hashVal = HashTim(Q2,W); break;
      //case kHarut: hashVal = HashHarut(hadE[kpip],hadE[kpim]); break;
      case kHarut: hashVal = evnum; break;
      case kSimple: hashVal = evnum; break;
    };

    hashMap.insert(std::pair<Int_t,Int_t>(hashVal,xi));
  };



  TString outdat = "compare.dat";
  gSystem->RedirectOutput(outdat,"w");
  gSystem->RedirectOutput(0);
  Bool_t extraCut;

  
  for(int i=0; i<ev->ENT; i++) {

    ev->GetEvent(i);

    if(ev->cutCrossCheck) {

      // hash tree's event
      switch(WHICH_XCHECK) {
        case kTim: hashVal = HashTim(ev->Q2,ev->W); break;
        //case kHarut: hashVal = HashHarut(ev->hadE[kpip],ev->hadE[kpim]); break;
        case kHarut: hashVal = ev->evnum; break;
        case kSimple: hashVal = ev->evnum; break;
      };

      // see if xtree has a matching hash value
      hashIter = hashMap.find(hashVal);
      if(hashIter!=hashMap.end()) {

        // set xtree variables to the matching event's
        xtree->GetEntry(hashIter->second);

        // extra requirement to improve event matching
        switch(WHICH_XCHECK) {
          case kTim: 
            extraCut = fabs(Q2 - ev->Q2) < 0.01;
            break;
          case kHarut: 
            //extraCut = 
              //fabs( hadTheta[kpip] - Tools::EtaToTheta(ev->hadEta[kpip]) ) < 0.005;
            extraCut = true;
            break;
          default: extraCut = true;
        };
        if(extraCut) {

          gSystem->RedirectOutput(outdat,"a");

          printf("EVENT#  xtree: %d  tree: %d\n",evnum,ev->evnum);
          //printf("HASH  xtree: %f  tree: %f\n",hashIter->first,hashVal);
          printf("%7s %7s %7s %7s\n","var","xtree","tree","|diff|");
          printf("%7s %7s %7s %7s\n","---","-----","----","------");

          switch(WHICH_XCHECK) {
            case kTim: 
              PrintCompare("Q2",Q2,ev->Q2);
              PrintCompare("W",W,ev->W);
              PrintCompare("x",x,ev->x);
              PrintCompare("y",y,ev->y);
              PrintCompare("Mh",Mh,ev->Mh);
              PrintCompare("xF",xF,ev->xF);
              PrintCompare("theta",TMath::Sin(theta),TMath::Sin(ev->theta));
              PrintCompare("PhiH",PhiH,ev->PhiH);
              PrintCompare("PhiR",PhiR,ev->PhiR);
              break;
            case kHarut:

              PrintCompare("pipE",hadE[kpip],ev->hadE[kpip]);
              PrintCompare("pipPt",hadPt[kpip],ev->hadPt[kpip]);
              //PrintCompare("pipPhi",hadPhi[kpip],ev->hadPhi[kpip]);
              //PrintCompare("pipTh",hadTheta[kpip],Tools::EtaToTheta(ev->hadEta[kpip]));

              PrintCompare("pimE",hadE[kpim],ev->hadE[kpim]);
              PrintCompare("pimPt",hadPt[kpim],ev->hadPt[kpim]);
              //PrintCompare("pimPhi",hadPhi[kpim],ev->hadPhi[kpim]);
              //PrintCompare("pimTh",hadTheta[kpim],Tools::EtaToTheta(ev->hadEta[kpim]));

              break;
            case kSimple:

              PrintCompare("pipE",hadE[kpip],ev->hadE[kpip]);
              PrintCompare("pipPt",hadPt[kpip],ev->hadPt[kpip]);

              PrintCompare("pimE",hadE[kpim],ev->hadE[kpim]);
              PrintCompare("pimPt",hadPt[kpim],ev->hadPt[kpim]);

              break;
          };
          printf("\n");

          gSystem->RedirectOutput(0);
        };

      };
    };
  };

  gROOT->ProcessLine(TString(".! cat "+outdat));

};


void PrintCompare(TString name, Float_t val, Float_t xval) {
  Float_t diff;

  // if it's an angle, ensure it's in proper range
  if(name.Contains("Phi")) {
    val = Tools::AdjAngleTwoPi(val);
    xval = Tools::AdjAngleTwoPi(xval);
    diff = Tools::AdjAngleTwoPi( fabs(val - xval) );
  } else {
    diff = fabs(val - xval);
  };

  // print comparison
  printf("%7s %7.2f %7.2f %7.2f\n",name.Data(),val,xval,diff);

};
