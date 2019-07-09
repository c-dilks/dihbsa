// for printing kinematics to cross check with others

R__LOAD_LIBRARY(../src/DihBsa)

#include "Constants.h"
#include "Tools.h"
#include "EventTree.h"
#include <map>
#include <math.h>


// hash function
Float_t Digitize(Float_t val) { 
  Int_t digits = 1; // number of decimal digits to keep
  return roundf(val*pow(10,digits));
};
Int_t HashTim(Float_t Q2_, Float_t W_) { 
  return 1000 * Digitize(Q2_) + Digitize(W_);
};
Int_t HashHarut(Float_t pipE_, Float_t pimE_) {
  return 1000 * Digitize(pipE_) + Digitize(pimE_);
};


const Int_t NF = 2;


// printing methods
void PrintCompare(TString name, Float_t val, Float_t xval);


void CrossChecker(TString indir="../outroot") {

  ///////////////////////////////
  enum xenum { kTim, kHarut, kSimpleC, kSimpleJava, kAnalysis };
  Int_t xcheck[NF] = { kSimpleC, kHarut };
  ///////////////////////////////


  // instantiate EventTree if using dihbsa analysis code (kAnalysis)
  EventTree * ev;
  if( xcheck[0]==kAnalysis || xcheck[1]==kAnalysis ) {
    ev = new EventTree( TString(indir+"/*.root"), EncodePairType(kPip,kPim) ); 
  };


  enum hadron_enum { iP, iM, nHad };
  TString hadN[nHad] = { "pip", "pim" };
  int h;



  ////////////////////////////
  // standard set of variables for comparison
  // -- IMPORTANT: make sure these are all reset in every iteration of the event loop !
  Int_t evnum[NF];
  Float_t Q2[NF], W[NF], x[NF], y[NF];
  Float_t Mh[NF], xF[NF], theta[NF], PhiH[NF], PhiR[NF];
  Float_t PhPt[NF];

  Float_t hadE[NF][nHad];
  Float_t hadPt[NF][nHad];
  Float_t hadTheta[NF][nHad];
  Float_t hadPhi[NF][nHad];
  Float_t hadPx[NF][nHad];
  Float_t hadPy[NF][nHad];
  Float_t hadPz[NF][nHad];
  ////////////////////////////

  Float_t valTheta;

  TFile * xfile[NF];
  TTree * xtree[NF];
  TString xstr;
  TString xtreeN[NF];
  int f;

  for(f=0; f<NF; f++) {

    if( f != kSimpleC ) {
      xtreeN[f] = Form("xtree%d",f);
      xtree[f] = new TTree(xtreeN[f],xtreeN[f]);
    };

    switch(xcheck[f]) {

    case kTim:

      gROOT->ProcessLine(".! python formatTimFile.py");
      xstr = "evnum/I:Q2/F:W/F:x/F:y/F:Mh/F:pT/F:xF/F:theta/F:PhiR/F:PhiH/F";
      printf("xtree[%d] branches: %s\n",f,xstr.Data());

      xtree[f]->ReadFile("xtree.dat",xstr);
      xtree[f]->SetBranchAddress("evnum",&evnum[f]);
      xtree[f]->SetBranchAddress("Q2",&Q2[f]);
      xtree[f]->SetBranchAddress("W",&W[f]);
      xtree[f]->SetBranchAddress("x",&x[f]);
      xtree[f]->SetBranchAddress("y",&y[f]);
      xtree[f]->SetBranchAddress("Mh",&Mh[f]);
      xtree[f]->SetBranchAddress("xF",&xF[f]);
      xtree[f]->SetBranchAddress("theta",&theta[f]);
      xtree[f]->SetBranchAddress("PhiH",&PhiH[f]);
      xtree[f]->SetBranchAddress("PhiR",&PhiR[f]);

      break;

    case kHarut:

      gROOT->ProcessLine(".! python formatHarutFile.py");
      xstr = "evnum/I";
      for(h=0; h<nHad; h++) {
        xstr += ":"+hadN[h]+"E/F";
        xstr += ":"+hadN[h]+"Theta/F";
        xstr += ":"+hadN[h]+"Phi/F";
        xstr += ":"+hadN[h]+"Pt/F";
      };
      xstr += ":PhPt/F";
      printf("xtree[%d] branches: %s\n",f,xstr.Data());
      xtree[f]->ReadFile("xtree.dat",xstr);

      xtree[f]->SetBranchAddress("evnum",&evnum[f]);
      xtree[f]->SetBranchAddress("PhPt",&PhPt[f]);
      for(h=0; h<nHad; h++) {
        xtree[f]->SetBranchAddress( TString(hadN[h]+"E"), &hadE[f][h] );
        xtree[f]->SetBranchAddress( TString(hadN[h]+"Pt"), &hadPt[f][h] );
        xtree[f]->SetBranchAddress( TString(hadN[h]+"Theta"), &hadTheta[f][h] );
        xtree[f]->SetBranchAddress( TString(hadN[h]+"Phi"), &hadPhi[f][h] );
      };

      break;

    case kSimpleC:

      xfile[f] = new TFile("../simpleTree.root","READ");
      xtree[f] = (TTree*) xfile[f]->Get("tree");

      xtree[f]->SetBranchAddress("evnum",&evnum[f]);
      for(h=0; h<nHad; h++) {
        xtree[f]->SetBranchAddress( TString(hadN[h]+"E"), &hadE[f][h] );
        xtree[f]->SetBranchAddress( TString(hadN[h]+"Pt"), &hadPt[f][h] );
        xtree[f]->SetBranchAddress( TString(hadN[h]+"Px"), &hadPx[f][h] );
        xtree[f]->SetBranchAddress( TString(hadN[h]+"Py"), &hadPy[f][h] );
        xtree[f]->SetBranchAddress( TString(hadN[h]+"Pz"), &hadPz[f][h] );
      };

      break;

    case kSimpleJava:
      
      xstr = "evnum/I";
      for(h=0; h<nHad; h++) {
        xstr += ":"+hadN[h]+"Px/F";
        xstr += ":"+hadN[h]+"Py/F";
        xstr += ":"+hadN[h]+"Pz/F";
        xstr += ":"+hadN[h]+"E/F";
        xstr += ":"+hadN[h]+"Pt/F";
      };
      printf("xtree[%d] branches: %s\n",f,xstr.Data());
      xtree[f]->ReadFile("../jsrc/javaOut.dat",xstr);

      xtree[f]->SetBranchAddress("evnum",&evnum[f]);
      for(h=0; h<nHad; h++) {
        xtree[f]->SetBranchAddress( TString(hadN[h]+"E"), &hadE[f][h] );
        xtree[f]->SetBranchAddress( TString(hadN[h]+"Pt"), &hadPt[f][h] );
        xtree[f]->SetBranchAddress( TString(hadN[h]+"Px"), &hadPx[f][h] );
        xtree[f]->SetBranchAddress( TString(hadN[h]+"Py"), &hadPy[f][h] );
        xtree[f]->SetBranchAddress( TString(hadN[h]+"Pz"), &hadPz[f][h] );
      };

      break;

    default: 
      fprintf(stderr,"ERROR: unknown xcheck[%d]=%d\n",f,xcheck[f]);
      exit(0);

    };
  };



  // build hash table for xtree[1]
  // -- later we will loop through xtree[0], searching for each event's hash value to
  //    this hash table in order to find a matching event in xtree[1]
  Int_t hashVal;
  std::map<Int_t,Int_t> hashMap; // event hash -> xtree[1] index
  std::map<Int_t,Int_t>::iterator hashIter;

  for(int xi=0; xi<xtree[1]->GetEntries(); xi++) {
    xtree[1]->GetEntry(xi);

    switch(xcheck[1]) {
      case kTim: hashVal = HashTim(Q2[1],W[1]); break;
      //case kHarut: hashVal = HashHarut(hadE[1][iP],hadE[1][iM]); break;
      default: hashVal = evnum[1];
    };

    hashMap.insert(std::pair<Int_t,Int_t>(hashVal,xi));
  };



  // define output file
  TString outdat = "compare.dat";
  gSystem->RedirectOutput(outdat,"w");
  gSystem->RedirectOutput(0);

  Bool_t extraCut,evCut;
  

  // loop through xtree[0]
  for(int i=0; i<xtree[0]->GetEntries(); i++) {

    // reset variables so that it's easy to filter out which
    // aren't associated to any branch
    for(f=0; f<NF; f++) {
      evnum[f] = -10000;
      Q2[f] = -10000;
      W[f] = -10000;
      x[f] = -10000;
      y[f] = -10000;
      Mh[f] = -10000;
      xF[f] = -10000;
      theta[f] = -10000;
      PhiH[f] = -10000;
      PhiR[f] = -10000;
      PhPt[f] = -10000;
      for(h=0; h<nHad; h++) {
        hadE[f][h] = -10000;
        hadPt[f][h] = -10000;
        hadTheta[f][h] = -10000;
        hadPhi[f][h] = -10000;
        hadPx[f][h] = -10000;
        hadPy[f][h] = -10000;
        hadPz[f][h] = -10000;
      };
    };

    // fill xtree[0] variables
    xtree[0]->GetEntry(i);

    //ev->GetEvent(i); // evtr loop
    //if(ev->cutCrossCheck) 


    evCut = true; // TODO
    if(evCut) {

      // hash xtree[0]'s event
      switch(xcheck[0]) {
        case kTim: hashVal = HashTim(Q2[0],W[0]); break;
        //case kHarut: hashVal = HashHarut(hadE[0][iP],hadE[0][iM]); break;
        default: hashVal = evnum[0];
      };


      // check if xtree[1] has a matching hash value
      hashIter = hashMap.find(hashVal);
      if(hashIter!=hashMap.end()) {

        // set xtree[1] variables to the matching event's
        xtree[1]->GetEntry(hashIter->second);

        // extra requirement to improve event matching
        /*
        switch(WHICH_XCHECK) {
          case kTim: 
            extraCut = fabs(Q2 - ev->Q2) < 0.01;
            break;
          case kHarut: 
            //extraCut = 
              //fabs( hadTheta[iP] - Tools::EtaToTheta(ev->hadEta[iP]) ) < 0.005;
            extraCut = true;
            break;
          default: extraCut = true;
        };
        */
        extraCut = true;
        if(extraCut) {

          // print comparisons
          // -- if a variable is not set (i.e., set to -10000), its comparison will not
          //    be printed
          gSystem->RedirectOutput(outdat,"a");

          printf("EVENT#");
          for(f=0; f<NF; f++) printf("  xtree%d: %d",f,evnum[f]); printf("\n");
          printf("HASH  xtree0: %d  xtree1: %d\n",hashVal,hashIter->first);
          printf("%8s %8s %8s %8s\n","var","xtree0","xtree1","|diff|");
          printf("%8s %8s %8s %8s\n","---","------","------","------");


          for(h=0; h<nHad; h++) {
            PrintCompare( TString(hadN[h]+"E"), hadE[0][h], hadE[1][h] );
            PrintCompare( TString(hadN[h]+"Pt"), hadPt[0][h], hadPt[1][h] );
            PrintCompare( TString(hadN[h]+"Px"), hadPx[0][h], hadPx[1][h] );
            PrintCompare( TString(hadN[h]+"Py"), hadPy[0][h], hadPy[1][h] );
            PrintCompare( TString(hadN[h]+"Pz"), hadPz[0][h], hadPz[1][h] );
            PrintCompare( TString(hadN[h]+"Theta"), hadTheta[0][h], hadTheta[1][h] );
            PrintCompare( TString(hadN[h]+"Phi"), hadPhi[0][h], hadPhi[1][h] );
          };

          PrintCompare( "Q2", Q2[0], Q2[1] );
          PrintCompare( "W", W[0], W[1] );
          PrintCompare( "x", x[0], x[1] );
          PrintCompare( "y", y[0], y[1] );
          PrintCompare( "Mh", Mh[0], Mh[1] );
          PrintCompare( "xF", xF[0], xF[1] );
          PrintCompare( "theta", TMath::Sin(theta[0]), TMath::Sin(theta[1]) );
          PrintCompare( "PhiH", PhiH[0], PhiH[1] );
          PrintCompare( "PhiR", PhiR[0], PhiR[1] );

          printf("\n");

          gSystem->RedirectOutput(0);
        };

      };
    };
  };

  gROOT->ProcessLine(TString(".! cat "+outdat));

};


void PrintCompare(TString name, Float_t val0, Float_t val1) {

  if(val0<-1000 || val1<-1000) return;

  Float_t diff;

  // if it's an angle, ensure it's in proper range
  if(name.Contains("Phi")) {
    val0 = Tools::AdjAngleTwoPi(val0);
    val1 = Tools::AdjAngleTwoPi(val1);
    diff = Tools::AdjAngleTwoPi( fabs(val0 - val1) );
  } else {
    diff = fabs(val0 - val1);
  };

  // print comparison
  printf("%8s %8.2f %8.2f %8.2f\n",name.Data(),val0,val1,diff);

};
