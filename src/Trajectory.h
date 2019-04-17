#ifndef Trajectory_
#define Trajectory_

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

// dihbsa
#include "Constants.h"



class Trajectory : public TObject
{
  public:
    Trajectory(Int_t particle_index);
    ~Trajectory();

    void SetMomentum(Float_t px, Float_t py, Float_t pz);
    TString Name() { return PartName(Idx); };
    TString Title() { return PartTitle(Idx); };
    Int_t PID() { return PartPID(Idx); };
    Float_t Mass() { return PartMass(Idx); };


    Int_t Idx;
    TLorentzVector * Vec;


    Bool_t debug;
    
  private:


  ClassDef(Trajectory,1);
};

#endif
