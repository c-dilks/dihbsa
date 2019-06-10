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
#include "Tools.h"




class EventTree : public TObject
{
  public:
    EventTree(TString filelist, Int_t whichPair_);
    ~EventTree();

    void GetEvent(Int_t i);
    Bool_t Valid();
    void PrintEvent();

    Bool_t debug;
    Long64_t ENT;


    ///////////////////////////
    //   BRANCHES
    ///////////////////////////
    // DIS kinematics
    Float_t W,Q2,Nu,x,y;

    // hadron kinematics
    Int_t hadIdx[2];
    Float_t hadE[2];
    Float_t hadP[2];
    Float_t hadPt[2];
    Float_t hadEta[2];
    Float_t hadPhi[2];

    // electron kinematics
    Float_t eleE;
    Float_t elePt;
    Float_t eleEta;
    Float_t elePhi;


    // dihadron kinematics
    Int_t particleCnt[nParticles];
    Int_t particleCntAll;
    Float_t Mh,Zpair,PhiH,Mmiss,xF,alpha;
    Float_t Z[2];
    Float_t Ph,PhPerp;
    Float_t PhEta,PhPhi;
    Float_t R,RPerp,RT;

    // event-level branches
    Int_t evnum,runnum;
    Int_t helicity;
    Float_t torus;
    Long64_t triggerBits;

    // PhiR 
    Float_t PhiR; // set to the preferred one
    Float_t PhiRq;
    Float_t PhiRp;
    Float_t PhiRp_r;
    Float_t PhiRp_g;

    Float_t b_PhiRq;
    Float_t b_PhiRp;
    Float_t b_PhiRp_r;
    Float_t b_PhiRp_g;

    // diphotons
    Float_t photE[2];
    Float_t photPt[2];
    Float_t photEta[2];
    Float_t photPhi[2];
    Float_t diphE;
    Float_t diphZ;
    Float_t diphPt;
    Float_t diphM;
    Float_t diphAlpha;
    Float_t diphEta;
    Float_t diphPhi;
    ///////////////////////////


    ///////////////////////////
    //   EventCuts
    ///////////////////////////
   Bool_t cutQ2,cutW,cutY,cutDIS;
   Bool_t cutDihadronKinematics;
   Bool_t cutDihadron;
   Bool_t cutDiphKinematics;
   Bool_t cutDiph;

   Bool_t useDiphBG;
    
  private:
    TChain * chain;
    Int_t whichHad[2];


  ClassDef(EventTree,1);
};

#endif
