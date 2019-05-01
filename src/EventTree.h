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
    Float_t Mh,Zpair,PhiH,PhiR,Mmiss,xF;
    Float_t Z[2];
    Float_t Ph,Pht;

    // event-level branches
    Int_t evnum,runnum;
    Int_t helicity;
    Float_t torus;
    Long64_t triggerBits;

    // PhiR tests
    Float_t PhiR_T_byKt;
    Float_t PhiR_T_byRej;
    Float_t PhiR_Perp;
    Float_t PhiR_byPh;
    Float_t PhiR_byPhad[2];
    Float_t PhiP1P2;

    Float_t b_PhiR_T_byKt;
    Float_t b_PhiR_T_byRej;
    Float_t b_PhiR_Perp;
    Float_t b_PhiR_byPh;
    Float_t b_PhiR_byPhad[2];
    Float_t b_PhiP1P2;
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
