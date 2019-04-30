#include "Dihadron.h"

ClassImp(Dihadron)

using namespace std;


Dihadron::Dihadron() {
  debug = true;
  useBreit = false;

  for(h=0; h<2; h++) vecHad[h] = TLorentzVector(0,0,0,0);
  vecPh = TLorentzVector(0,0,0,0);
  vecR = TLorentzVector(0,0,0,0);

  printf("Dihadron instantiated\n");
};


void Dihadron::SetEvent(
  Trajectory * trajPlus, 
  Trajectory * trajMinus,
  DIS * disEvent
) {
  hadron[hP] = trajPlus;
  hadron[hM] = trajMinus;
  disEv = disEvent;


  // get disEv vectors
  disVecBeam = disEv->vecBeam;
  disVecTarget = disEv->vecTarget;
  disVecElectron = disEv->vecElectron;
  disVecW = disEv->vecW;
  disVecQ = disEv->vecQ;
  if(useBreit) {
    disVecBeam.Boost(disEv->BreitBoost);
    disVecTarget.Boost(disEv->BreitBoost);
    disVecElectron.Boost(disEv->BreitBoost);
    disVecW.Boost(disEv->BreitBoost);
    disVecQ.Boost(disEv->BreitBoost);
  };


  // compute 4-momenta Ph and R
  for(h=0; h<2; h++) {
    vecHad[h] = hadron[h]->Vec;
    if(useBreit) vecHad[h].Boost(disEv->BreitBoost);
  };
  vecPh = vecHad[hP] + vecHad[hM];
  vecR = 0.5 * ( vecHad[hP] - vecHad[hM] );

  // compute z
  for(h=0; h<2; h++) z[h] = vecHad[h].E() / disVecQ.E();
  zpair = vecPh.E() / disVecQ.E();

  // compute invariant mass
  Mh = vecPh.M();

  // compute missing mass
  vecMmiss = disVecW - vecPh;
  Mmiss = vecMmiss.M();

  // compute xF
  bvecPh = vecPh;
  if(!useBreit) bvecPh.Boost(disEv->BreitBoost); // if useBreit==true, bvecPh is already
                                                 // in the Breit frame
  xF = bvecPh.E() / disEv->W;


  // compute angles
  ComputeAngles();
};


void Dihadron::ComputeAngles() {
  // sign convention:
  //   P_1 = pi+ momentum
  //   P_2 = pi- momentum
  //  
  // regarding transverse components, there are two 
  // frames to consider (see arXiv:1707.04999):
  // -- perp-frame: "transverse" plane is normal to fragmenting
  //                quark, i.e., to q
  // -- T-frame: "transverse" plane is normal to Ph


  // get 3-momenta from 4-momenta
  pQ = disVecQ.Vect();
  pL = disVecElectron.Vect();
  pPh = vecPh.Vect();
  pR = vecR.Vect();
  for(h=0; h<2; h++) pHad[h] = vecHad[h].Vect();


  // perp-frame hadron 3-momenta
  for(h=0; h<2; h++) pHad_Perp[h] = Reject(pHad[h],pQ);
  pPh_Perp = Reject(pPh,pQ);


  // compute transverse components of R in different frames
  // (and with different equations)
  // -- in T-frame, via kT relation (following 1707.04999)
  pR_T_byKt = 1 / ( z[hP] + z[hM] ) *  
              ( z[hM] * pHad_Perp[hP]  -  z[hP] * pHad_Perp[hM] );
  // -- in T-frame, via rejection
  pR_T_byRej = Reject(pR,pPh);
  // -- in perp-frame, via rejection
  pR_Perp= Reject(pR,pQ);


  // momentum magnitudes
  PhMag = pPh.Mag(); // dihadron total momentum Ph
  PhtMag = pPh_Perp.Mag(); // trans. comp. of Ph (used for G1perp)
  RMag = vecR.Mag(); // dihadron relative momentum R
  RtMag = pR_T_byKt.Mag(); // trans. comp. of R


  // compute PhiH angle
  PhiH = PlaneAngle(pQ,pL,pQ,pPh);


  // compute PhiR angle (all the ways)

  // -- HERMES 0803.2367 angle, but used Matevosyan et al 1707.04999
  //    to obtain R_T vector
  PhiR_T_byKt = PlaneAngle(pQ,pL,pQ,pR_T_byKt);

  // -- HERMES 0803.2367 angle
  PhiR_T_byRej = PlaneAngle(pQ,pL,pQ,pR_T_byRej);

  // -- COMPASS 1702.07317, but used vector rejection to get R_perp
  PhiR_Perp = PlaneAngle(pQ,pL,pQ,pR_Perp);

  // -- alternative tests
  PhiR_byPh = PlaneAngle(pQ,pL,pQ,pPh);
  PhiP1P2 = PlaneAngle(pQ,pL,pHad[hM],pHad[hP]);
  for(int h=0; h<2; h++) 
    PhiR_byPhad[h] = PlaneAngle(pQ,pL,pQ,pHad[h]);

  // finally set the "preferred" PhiR angle
  PhiR = PhiR_T_byKt;
};


// compute angle between planes given by
// vA x vB and vC x vD
// the angle is with respect to vA x vB
Float_t Dihadron::PlaneAngle(
  TVector3 vA, TVector3 vB,
  TVector3 vC, TVector3 vD
) {
  crossAB = vA.Cross(vB); // AxB
  crossCD = vC.Cross(vD); // CxD

  sgn = crossAB.Dot(vD); // (AxB).D

  if(fabs(sgn)<0.0001) {
    //fprintf(stderr,"WARNING: Dihadron::PlaneAngle (AxB).D == 0\n");
    return -10000;
  };

  sgn /= fabs(sgn); // sign of (AxB).D

  numer = crossAB.Dot(crossCD); // (AxB).(CxD)
  denom = crossAB.Mag() * crossCD.Mag(); // |AxB|*|CxD|

  if(fabs(denom)<0.0001) {
    //fprintf(stderr,"WARNING: Dihadron::PlaneAngle |AxB|*|CxD| == 0\n");
    return -10000;
  };

  return sgn * TMath::ACos(numer/denom);
};


// vector rejection: 
// returns vA projected onto plane transverse to vB
TVector3 Dihadron::Reject(TVector3 vA, TVector3 vB) {

  if(fabs(vB.Dot(vB))<0.0001) {
    //fprintf(stderr,"WARNING: Dihadron::Reject to null vector\n");
    return TVector3(0,0,0);
  };

  return vA - Project(vA,vB);

};

// vector projection:
// returns vA projected onto vB
TVector3 Dihadron::Project(TVector3 vA, TVector3 vB) {

  if(fabs(vB.Dot(vB))<0.0001) {
    //fprintf(stderr,"WARNING: Dihadron::Project to null vector\n");
    return TVector3(0,0,0);
  };

  proj = vA.Dot(vB) / ( vB.Dot(vB) );
  return proj * vB;

};

  

Dihadron::~Dihadron() {
};

