#ifndef EventTree_
#define EventTree_

#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <map>
#include <vector>

// ROOT
#include "TSystem.h"
#include "TObject.h"
#include "TTree.h"
#include "TChain.h"
#include "TFile.h"
#include "TLorentzVector.h"
#include "TString.h"
#include "TMath.h"

// dihbsa
#include "Constants.h"
#include "Trajectory.h"
#include "DIS.h"




class EventTree : public TObject
{
  public:
    EventTree(TString filelist);
    ~EventTree();

    void GetEvent(Int_t i);
    void Print();

    Bool_t debug;
    Int_t ENT;


    ///////////////////////////
    //   BRANCHES
    ///////////////////////////
    // DIS kinematics
    Float_t W,Q2,Nu,x,y;

    // hadron kinematics
    Float_t hadE[2]; // [enum plus_minus (0=+, 1=-)]
    Float_t hadP[2];
    Float_t hadPt[2];
    Float_t hadEta[2];
    Float_t hadPhi[2];

    // dihadron kinematics
    Float_t Mh,Zpair,PhiH,PhiR,Mmiss,xF,alpha;
    Float_t Z[2];
    Float_t Ph,PhPerp;
    Float_t R,RPerp,RT;

    // event-level branches
    Int_t evnum,runnum;
    Int_t helicity;
    Float_t torus;
    Long64_t triggerBits;

    // PhiR tests
    Float_t PhiRq;
    Float_t PhiRp;
    Float_t PhiRp_r;

    Float_t b_PhiRq;
    Float_t b_PhiRp;
    Float_t b_PhiRp_r;
    ///////////////////////////


    ///////////////////////////
    //   EventCuts
    ///////////////////////////
   Bool_t cutQ2,cutW,cutY;
   Bool_t cutDihadron;
    
  private:
    TChain * chain;


  ClassDef(EventTree,1);
};

#endif
