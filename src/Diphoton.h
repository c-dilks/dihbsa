#ifndef Diphoton_
#define Diphoton_

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



class Diphoton : public TObject
{
  public:
    Diphoton();
    ~Diphoton();

    Bool_t debug;

    void SetEvent(Trajectory * traj1, Trajectory * traj2);
    void ResetVars();

    Trajectory * photon[2]; // photon trajectories
    Trajectory * Traj; // diphoton trajectory


    // diphoton variables
    Float_t E; // energy of the pair
    Float_t Ephot[2]; // energy of each photon
    Float_t Z; // energy sharing (E1-E2)/E
    Float_t Pt; // transverse momentum
    Float_t Mgg; // invariant mass
    Float_t Alpha; // diphoton opening angle
    Float_t Eta; // pseudorapidity
    Float_t Phi; // azimuth

    // booleans
    Bool_t validDiphoton;



  private:

    TLorentzVector vecDiphoton;
    TVector3 momPhoton[2];



  ClassDef(Diphoton,1);
};

#endif
