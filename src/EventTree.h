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
#include "TString.h"
#include "TMath.h"
#include "TRandom.h"

// dihbsa
#include "Constants.h"
#include "Trajectory.h"
#include "DIS.h"
#include "Dihadron.h"
#include "Tools.h"




class EventTree : public TObject
{
  public:
    EventTree(TString filelist, Int_t whichPair_);
    ~EventTree();

    void GetEvent(Int_t i);
    Int_t SpinState();
    Bool_t Valid();
    void PrintEvent();

    Dihadron * GetDihadronObj();
    DIS * GetDISObj();

    Bool_t BuildEvnumMap();
    Bool_t FindEvent(Int_t evnum_);

    Bool_t debug;
    Long64_t ENT;


    ///////////////////////////
    //   BRANCHES
    ///////////////////////////
    // DIS kinematics
    Float_t W,Q2,Nu,x,y;

    // hadron kinematics
    Int_t pairType;
    Int_t hadIdx[2];
    Float_t hadE[2];
    Float_t hadP[2];
    Float_t hadPt[2];
    Float_t hadPtq[2];
    Float_t hadEta[2];
    Float_t hadPhi[2];
    Float_t hadXF[2];

    // electron kinematics
    Float_t eleE;
    Float_t eleP;
    Float_t elePt;
    Float_t eleEta;
    Float_t elePhi;


    // dihadron kinematics
    Int_t particleCnt[nParticles];
    Int_t particleCntAll;
    Float_t Mh,Zpair,PhiH,Mmiss,xF,alpha;
    Float_t Z[2];
    Float_t zeta;
    Float_t theta;

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

    Float_t PhiHR; // PhiH-PhiR

    /*
    Float_t b_PhiRq;
    Float_t b_PhiRp;
    Float_t b_PhiRp_r;
    Float_t b_PhiRp_g;
    */

    // diphotons
    Int_t diphCnt;
    Float_t diphPhotE[2][2]; // [diphCnt] [photon 0 or 1]
    Float_t diphPhotPt[2][2];
    Float_t diphPhotEta[2][2];
    Float_t diphPhotPhi[2][2];
    Float_t diphE[2]; // [diphCnt]
    Float_t diphZ[2];
    Float_t diphPt[2];
    Float_t diphM[2];
    Float_t diphAlpha[2];
    Float_t diphEta[2];
    Float_t diphPhi[2];
    Float_t angEle[2][2];
    ///////////////////////////


    ///////////////////////////
    //   EventCuts
    ///////////////////////////
   Bool_t cutQ2,cutW,cutY,cutDIS;
   Bool_t cutDihadronKinematics;
   Bool_t cutDihadron;
   Bool_t cutDiphKinematics[2];
   Bool_t cutDiph[2];
   Bool_t cutCrossCheck;

   Bool_t useDiphBG;
    
  private:
    TChain * chain;
    Int_t whichHad[2];

    TRandom * RNG;

    Dihadron * objDihadron;
    DIS * objDIS;
    Trajectory * trEle;
    Trajectory * trHad[2];

    TLorentzVector photMom[2];
    TLorentzVector hadMom[2];
    TLorentzVector eleMom;

    std::map<Int_t,Int_t> evnumMap;
    bool inserted;
    std::map<Int_t,Int_t>::iterator foundIt;


  ClassDef(EventTree,1);
};

#endif
