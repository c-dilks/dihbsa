#ifndef Dihadron_
#define Dihadron_

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
#include "TFile.h"
#include "TLorentzVector.h"
#include "TString.h"
#include "TMath.h"

// dihbsa
#include "Constants.h"
#include "Trajectory.h"
#include "DIS.h"




class Dihadron : public TObject
{
  public:
    Dihadron();
    ~Dihadron();

    Bool_t debug;

    void SetEvent(
      Trajectory * trajPlus, Trajectory * trajMinus, DIS * disEv);
    void ComputeAngles();
    Float_t PlaneAngle(TVector3 vA, TVector3 vB,
                       TVector3 vC, TVector3 vD);
    TVector3 Reject(TVector3 vA, TVector3 vB);
    TVector3 Project(TVector3 vA, TVector3 vB);

    Trajectory * hadron[2];
    DIS * disEv;

    TLorentzVector vecHad[2]; // hadron momentum
    TLorentzVector vecPh; // dihadron total momentum
    TLorentzVector vecR; // dihadron relative momentum
    TLorentzVector vecMmiss; // used to compute missing mass

    TLorentzVector bvecPh; // breit frame's dihadron total momentum


    Float_t PhMag; // dihadron total momentum
    Float_t PhtMag; // transverse component of dihadron total momentum
    Float_t RMag; // dihadron relative momentum
    Float_t RtMag; // transverse componet of relative momentum
    Float_t phiR; // angle[ reaction_plane, R^q etc. (see below) ]
    Float_t phiH; // angle[ reaction_plane, Ph^q ]
    Float_t z[2]; // fraction of energy of fragmenting parton
                  // carried by the hadron
    Float_t zpair; // fraction of energy of fragmenting parton
                   // carried by the hadron pair
    Float_t Mh; // dihadron invariant mass
    Float_t Mmiss; // missing mass
    Float_t xF; // feynman-x


    // phiR angle is defined a couple different ways; these
    // variables are for those and a couple alternative tests
    Float_t phiR_T_byKt;
    Float_t phiR_T_byRej;
    Float_t phiR_Perp;
    Float_t phiR_byPh;
    Float_t phiR_byPhad[2];
    Float_t phiP1P2;

    
  private:
    int h;

    TVector3 pQ,pL,pPh,pR;
    TVector3 pHad[2];
    TVector3 pHad_Perp[2];
    TVector3 pPh_Perp;

    TVector3 pR_T_byKt;
    TVector3 pR_T_byRej;
    TVector3 pR_Perp;

    Float_t proj;
    TVector3 crossAB,crossCD;
    Float_t sgn,numer,denom;



  ClassDef(Dihadron,1);
};

#endif
