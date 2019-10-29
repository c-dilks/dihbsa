// plots asymmetry modulation versus various other kinematics, to get an overall sense
// of the acceptance from the data

R__LOAD_LIBRARY(DihBsa)

#include "Constants.h"
#include "EventTree.h"
#include "Tools.h"
#include "TString.h"
#include "TTree.h"

Int_t nTotal;
void PrintFraction(Int_t n, TString s);

void PrintNumEvents(TString fname=
"outroot.MC.rec/clasdispr.00.e10.600.emn0.75tmn.09.xs81.61nb.dis.0410.dat.hipo.root"
) {

  Int_t pairType = EncodePairType(kPip,kPim); // make sure ordering obeys dihHadIdx...
  EventTree * ev = new EventTree(fname,pairType);

  nTotal = 0;
  Int_t nValid = 0;
  Int_t nDIS = 0;
  Int_t nDihadronKin = 0;
  Int_t nVertex = 0;
  Int_t nFiducial = 0;



  for(int i=0; i<ev->ENT; i++) {
    ev->GetEvent(i);

    // check dihadron cuts
    if(pairType==ev->pairType) {
      nTotal++;
      if(ev->cutDIS) nDIS++;
      if(ev->cutDihadronKinematics) nDihadronKin++;
      if(ev->cutVertex) nVertex++;
      if(ev->cutFiducial) nFiducial++;
      if(ev->Valid()) nValid++;
    };
  };

  printf("\n\n");
  PrintFraction(nDIS,"DIS");
  PrintFraction(nDihadronKin,"dihadron kinematic");
  PrintFraction(nVertex,"vertex");
  PrintFraction(nFiducial,"e- fiducial cuts");
  PrintFraction(nValid,"full");
  printf("\n\n");

};

void PrintFraction(Int_t n, TString s) {
  printf("%d / %d dihadrons pass %s cuts (%.3f%%)\n",
    n,nTotal,s.Data(),Float_t(100*n)/nTotal);
};
