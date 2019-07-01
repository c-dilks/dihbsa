// for printing kinematics to cross check with Tim

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
Int_t Hash(Float_t Q2_,Float_t W_) { 
  return 1000 * Digitize(Q2_) + Digitize(W_);
};


// printing methods
void PrintCompare(TString name, Float_t val, Float_t xval);


void CrossCheckTim(TString indir="../outroot.crosscheck") {

  Int_t whichPair = EncodePairType(kPip,kPim);
  EventTree * ev = new EventTree(TString(indir+"/*.root"),whichPair); 


  TTree * xtree = new TTree();

  Int_t evnum;
  Float_t Q2,W,x,y,Mh,pT,xF,theta,PhiH,PhiR;
  Float_t elePol,XelePol;
  Float_t eleEta,XeleEta;

  // tim's tree
  gROOT->ProcessLine(".! python formatTimFile.py");
  xtree->ReadFile("xtree.dat","evnum/I:Q2/F:W/F:x/F:y/F:Mh/F:pT/F:xF/F:theta/F:PhiR/F:PhiH/F");
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
  //xtree->Scan("*");

  TString outdat = "compare.dat";
  gSystem->RedirectOutput(outdat,"w");
  gSystem->RedirectOutput(0);


  Float_t hashVal;
  std::map<Float_t,Int_t> hashMap;
  std::map<Float_t,Int_t>::iterator hashIter;
  for(int xi=0; xi<xtree->GetEntries(); xi++) {
    xtree->GetEntry(xi);
    hashVal = Hash(Q2,W);
    hashMap.insert(std::pair<Int_t,Int_t>(hashVal,xi));
  };


  
  for(int i=0; i<ev->ENT; i++) {

    ev->GetEvent(i);

    if(ev->pairType == whichPair) {

      hashVal = Hash(ev->Q2,ev->W);

      // see if xtree has a matching hash value
      hashIter = hashMap.find(hashVal);
      if(hashIter!=hashMap.end()) {

        // set xtree variables to the matching event's
        xtree->GetEntry(hashIter->second);

        // extra requirement to improve event matching
        if(fabs(Q2 - ev->Q2) < 0.01) {

          gSystem->RedirectOutput(outdat,"a");

          printf("EVENT#  xtree: %d  tree: %d\n",evnum,ev->evnum);
          //printf("HASH  xtree: %f  tree: %f\n",hashIter->first,hashVal);
          printf("%7s %7s %7s %7s\n","var","xtree","tree","|diff|");
          printf("%7s %7s %7s %7s\n","---","-----","----","------");
          PrintCompare("Q2",Q2,ev->Q2);
          PrintCompare("W",W,ev->W);
          PrintCompare("x",x,ev->x);
          PrintCompare("y",y,ev->y);
          PrintCompare("Mh",Mh,ev->Mh);
          //PrintCompare("pT",pT,ev->pT);
          PrintCompare("xF",xF,ev->xF);

          PrintCompare("theta",TMath::Sin(theta),TMath::Sin(ev->theta));

          PrintCompare("PhiH",PhiH,ev->PhiH);
          PrintCompare("PhiR",PhiR,ev->PhiR);

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
  if(name=="PhiH" || name=="PhiR") {
    val = Tools::AdjAngleTwoPi(val);
    xval = Tools::AdjAngleTwoPi(xval);
    diff = Tools::AdjAngleTwoPi( fabs(val - xval) );
  } else {
    diff = fabs(val - xval);
  };

  printf("%7s %7.2f %7.2f %7.2f\n",name.Data(),val,xval,diff);

};
