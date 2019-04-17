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

// dihbsa
#include "Constants.h"
#include "Trajectory.h"



class Dihadron : public TObject
{
  public:
    Dihadron();
    ~Dihadron();

    Bool_t debug;

    void SetHadrons(Trajectory * trajPlus, Trajectory * trajMinus);
    Float_t Mh();

    enum plusminus {hP,hM};
    Trajectory * hadron[2];
    TLorentzVector vecHad[2];
    TLorentzVector vecPh;
    TLorentzVector vecR;
    
  private:
    int h;


  ClassDef(Dihadron,1);
};

#endif
