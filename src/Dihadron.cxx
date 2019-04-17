#include "Dihadron.h"

ClassImp(Dihadron)

using namespace std;


Dihadron::Dihadron() {
  debug = true;

  for(h=0; h<2; h++) vecHad[h] = TLorentzVector(0,0,0,0);
  vecPh = TLorentzVector(0,0,0,0);
  vecR = TLorentzVector(0,0,0,0);

  printf("Dihadron instantiated\n");
};


void Dihadron::SetEvent(
  Trajectory * trajPlus, 
  Trajectory * trajMinus,
  DIS * disEv
) {
  hadron[hP] = trajPlus;
  hadron[hM] = trajMinus;
  disKin = disEv;

  // compute momenta
  for(h=0; h<2; h++) vecHad[h] = *(hadron[h]->Vec);
  vecPh = vecHad[hP] + vecHad[hM];
  vecR = 0.5 * ( vecHad[hP] - vecHad[hM] );

  // compute z
  for(h=0; h<2; h++) z[h] = hadron[h].Vec.E() / disKin->Nu;
  zpair = vecPh.E() / disKin->Nu;


  // compute angles
  ComputeAngles();
};


void Dihadron::ComputeAngles() {
  // using arxiv:1702.07317, with sign convention
  //   P_1 = pi+ momentum
  //   P_2 = pi- momentum
  //  

  // -- 3-momenta
  pQ = (disEv->vecQ).Vect();
  pL = (disEv->vecElectron).Vect();
  pPh = vecPh.Vect();
  pR = vecR.Vect();
  for(h=0; h<2; h++) {
    pHad[h] = vecHad[h].Vect();
    pHadPerp[h] = TVector3(pHad[h].X(),pHad[h].Y(),0);
  };

  // -- Rperp from equation 5
  pRperp = 1 / ( z[hP] + z[hM] ) * 
           ( z[hM] * pHadPerp[hP]  -  z[hP] * pHadPerp[hM] );

  // -- cross products
  crossQL = pQ.Cross(pL); // Q x l
  crossQPh = pQ.Cross(pPh); // Q x P_h
  crossQRperp = pQ.Cross(pRperp); // Q x Rperp

  // -- sign coefficients from equations 3 & 4
  sgnH = crossQL.Dot(pPh);
  sgnH /= fabs(sgnH); // sign of (Qxl).Ph
  sgnR = crossQL.Dot(pRperp);
  sgnR /= fabs(sgnR); // sign of (Qxl).Rperp

  // -- numerator and denominator of arcsin argument from eq 3 & 4
  numerH = crossQL.Dot(crossQPh); // (Qxl).(QxPh)
  denomH = crossQL.Mag() * crossQPh.Mag();

  numerR = crossQL.Dot(crossQRperp); // (Qxl).(QxRperp)
  denomR = crossQL.Mag() * crossQRperp.Mag();

  // -- finally compute the angles
  phiH = sgnH * TMath::ACos( numerH / denomH );
  phiR = sgnR * TMath::ACos( numerR / denomR );

};


Float_t Dihadron::Mh() {
  return vecPh.M();
};




Dihadron::~Dihadron() {
};

