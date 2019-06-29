// for printing kinematics to cross check with Tim

R__LOAD_LIBRARY(../src/DihBsa)

#include "Constants.h"
#include "EventTree.h"

void Compare(TString name, Float_t val, Float_t xval);

void CrossCheckTim(TString indir="../outroot.crosscheck") {

  Int_t whichPair = EncodePairType(kPip,kPim);
  EventTree * ev = new EventTree(TString(indir+"/*.root"),whichPair); 


  TTree * xtree = new TTree();

  // tim's tree
  xtree->ReadFile("xcheck.dat","evnum/I:Q2/F:W/F:x/F:y/F:Mh/F:pT/F:xF/F:theta/F:PhiR/F:PhiH/F");
  gROOT->ProcessLine(".! python formatTimFile.py");
  Int_t evnum;
  Float_t Q2,W,x,y,Mh,pT,xF,theta,PhiH,PhiR;
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

  



  
  for(int i=0; i<ev->ENT; i++) {

    ev->GetEvent(i);
    //ev->PrintEvent();

    if(ev->pairType == whichPair) {

      for(int xi=0; xi<xtree->GetEntries(); xi++) {
        xtree->GetEntry(xi);
        //printf("evnum=%d\n",evnum);

        if(ev->evnum +1 == evnum) {
          printf("EVENT# xtree=%d tree=%d\n",evnum,ev->evnum);
          printf("var\txtree\ttree\tdiff\n");
          printf("---\t-----\t----\t----\n");
          Compare("Q2",Q2,ev->Q2);
          Compare("W",W,ev->W);
          Compare("x",x,ev->x);
          Compare("y",y,ev->y);
          Compare("Mh",Mh,ev->Mh);
          //Compare("pT",pT,ev->pT);
          Compare("xF",xF,ev->xF);
          Compare("theta",theta,ev->theta);
          Compare("PhiH",PhiH,ev->PhiH);
          Compare("PhiR",PhiR,ev->PhiR);
          printf("---------\n");
        };
      };
    };
  };
};


void Compare(TString name, Float_t val, Float_t xval) {
  printf("%s:\t%.2f\t%.2f\t%.2f\n",name.Data(),val,xval,val-xval);
};


       

