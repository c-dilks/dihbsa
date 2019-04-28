#ifndef Asymmetry_
#define Asymmetry_

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
#include "TH1.h"

// dihbsa
#include "Constants.h"
#include "Trajectory.h"
#include "DIS.h"




class Asymmetry : public TObject
{
  public:
    Asymmetry();
    ~Asymmetry();
    void AddBinBound(Int_t iv, Float_t newBound);
    void PrintBinBounds();

    enum ivEnum { vM, vX, vZ, nIV }; // Independent Variables (IV)
    enum spinEnum { sU, sD, nSpin }; // spin state
    static const Int_t nBinsMax = 20;

    // 1-dim binning:
    // distributions for each bin for a specific IV; 
    // [spin] [IV] [bin]
    // -- [IV] is that specific IV; the other 2 IVs are integrated 
    // -- [bin] is the bin in that IV
    TH1D * pDist1[nSpin][nIV][nBinsMax]; // filled with, e.g., Sin(phiR) for each spin
    TH1D * wDist1[nIV][nBinsMax]; // filled with appropriate IV (finely binned)

    // 2-dim binning:
    // distributions for each pair of bins in 2 IVs; the third
    // IV not considered is integrated over
    // [spin] [IV] [bin1] [bin2]
    // -- [IV] is the IV that is integrated over
    // -- [bin1] is the first IV's bin
    // -- [bin2] is the second IV's bin
    // --> by convention, bin1's IV is in alphabetical order w.r.t.
    //     bin2's IV; e.g., [vX][0][1] is the 0th vM bin and 1st vZ bin
    TH1D * pDist2[nSpin][nIV][nBinsMax][nBinsMax];
    TH1D * wDist2[nIV][nBinsMax][nBinsMax];

    // 3-dim binning:
    // distributions for each triple of bins; since no IVs are integated
    // over, we only need one distribution; the three bins are in alphabetical
    // order
    // [vM bin] [vX bin] [vZ bin]
    TH1D * pDist3[nSpin][nBinsMax][nBinsMax][nBinsMax];
    TH1D * wDist3[nBinsMax][nBinsMax][nBinsMax];


    // bin boundaries
    Int_t nBins[nIV]; // the number of bins
    Float_t bound[nIV][nBinsMax]; // [IV] [bin#]; this is the lower bound
                                  // for bin "bin#"; to get the upper bound
                                  // of the highest bin, use bin#+1

    
  private:
    TString IVname[nIV];
    TString IVtitle[nIV];


  ClassDef(Asymmetry,1);
};

#endif
