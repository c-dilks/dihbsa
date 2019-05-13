#ifndef DIS_
#define DIS_

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

// dihbsa
#include "Constants.h"



class DIS : public TObject
{
  public:
    DIS();
    ~DIS();

    void SetBeamEn(Float_t newBeamEn);
    void SetElectron(Float_t px, Float_t py, Float_t pz);
    void Analyse();
    void PrintEvent();

    void ComputeBreitKinematics();
    void ComputeComKinematics();

    void BreitPrint();
    void ComPrint();

    void ResetVars();


    Float_t BeamEn;
    Float_t W,Q2,Nu,x,y;
    TVector3 BreitBoost;
    TVector3 ComBoost;

    Bool_t debug,speedup;
    
    // lab frame 4-vectors
    TLorentzVector vecBeam;
    TLorentzVector vecTarget;
    TLorentzVector vecElectron;
    TLorentzVector vecW;
    TLorentzVector vecQ;


    // breit frame 4-vectors 
    // (computed only if ComputeBreitKinematics() is called; 
    //  in practice, it's faster to externally boost only the ones you need)
    TLorentzVector breitBeam;
    TLorentzVector breitTarget;
    TLorentzVector breitElectron;
    TLorentzVector breitW;
    TLorentzVector breitQ;

    // target-virualPhoton center-of-mass frame 4-vectors 
    // (computed only if ComputeComKinematics() is called; 
    //  in practice, it's faster to externally boost only the ones you need)
    TLorentzVector comBeam;
    TLorentzVector comTarget;
    TLorentzVector comElectron;
    TLorentzVector comW;
    TLorentzVector comQ;

  private:
    TLorentzVector boostvecBreit;
    TLorentzVector boostvecCom;

  ClassDef(DIS,1);
};

#endif
