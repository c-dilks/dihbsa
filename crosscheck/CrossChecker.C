// for printing kinematics to cross check with others

R__LOAD_LIBRARY(../libDihBsa)

#include "Constants.h"
#include "Tools.h"
#include "EventTree.h"
#include <map>
#include <math.h>

const Int_t NF = 2; // number of files


///////////////////////////////
//
enum xenum { kTim, kHarut, kHarutOS, kOrlando, 
             kSimpleC, kSimpleJava, kAnalysis, kAnalysisLund };

Int_t xcheck[NF] = { kAnalysis, kTim };
//Int_t xcheck[NF] = { kAnalysisLund, kOrlando };
//Int_t xcheck[NF] = { kAnalysisLund, kHarut };
//Int_t xcheck[NF] = { kAnalysisLund, kHarutOS };

//Int_t xcheck[NF] = { kHarut, kOrlando };

//
///////////////////////////////


EventTree * ev[NF];
TFile * xfile[NF];
TTree * xtree[NF];
TString xstr;
Int_t ENT[NF];
enum hadron_enum { iP, iM, nHad };
TString hadN[nHad] = { "pip", "pim" };


////////////////////////////
// standard set of variables for comparison
// -- IMPORTANT: make sure these are all reset in every iteration of the event loop !
Int_t evnum[NF];
Float_t Q2[NF], W[NF], x[NF], y[NF];
Float_t Mh[NF], xF[NF], theta[NF], PhiH[NF], PhiR[NF];
Float_t PhPerp[NF];
Float_t eleP[NF];

Float_t hadE[NF][nHad];
Float_t hadP[NF][nHad];
Float_t hadPt[NF][nHad];
Float_t hadPtq[NF][nHad];
Float_t hadTheta[NF][nHad];
Float_t hadPhi[NF][nHad];
Float_t hadPx[NF][nHad];
Float_t hadPy[NF][nHad];
Float_t hadPz[NF][nHad];
Float_t hadZ[NF][nHad];
Float_t hadXF[NF][nHad];
////////////////////////////



void PrintCompare(TString name, Float_t val, Float_t xval);
void GetTreeVars(Int_t ff, Int_t ii);
void GetEventTreeVars(Int_t ff, Int_t ii);
void ReadOrlandoFormat(TString datname, Int_t ff);


void CrossChecker() {


  int h,f;

  // instantiate EventTree if using dihbsa analysis.exe
   for(f=0; f<NF; f++) {
    switch(xcheck[f]) {
      case kAnalysis:
        ev[f] = new EventTree( "../outrootCC/rec.root", EncodePairType(kPip,kPim) ); 
        break;
      case kAnalysisLund:
        ev[f] = new EventTree( "../outrootCC/lund.root", EncodePairType(kPip,kPim) ); 
        break;
    };
  };

  TString xtreeN[NF];

  char sep[20]; strcpy(sep,"--------");
  TString compareTitle[NF];
  TString compareName[NF];
  for(f=0; f<NF; f++) {

    if( f != kSimpleC ) {
      xtreeN[f] = Form("xtree%d",f);
      xtree[f] = new TTree(xtreeN[f],xtreeN[f]);
    };

    switch(xcheck[f]) {

      case kTim:
        compareTitle[f] = "Timothy -- hayward_cross_check_os_format.txt";
        compareName[f] = "timothy";
        gROOT->ProcessLine(
          ".! tail -n +2 xfiles/hayward_cross_check_os_format.txt > xtreeTim.dat"
        ); // (strip header)
        xstr = "evnum/I:Q2/F:x/F:W/F:y/F:eleP/F";
        xstr += ":Mh/F:xF/F:PhPerp/F:theta/F";
        xstr += ":PhiR/F:PhiH/F:pipP/F:pimP/F";

        /*
        // old file
        gROOT->ProcessLine(".! python formatTimFile.py");
        xstr = "evnum/I:eleE/F:pipE/F:pimE/F";
        xstr += ":Q2/F:W/F:x/F:y/F:Mh/F:pT/F:xF/F:theta/F:PhiR/F:PhiH/F";
        */

        printf("xtree[%d] branches: %s\n",f,xstr.Data());
        xtree[f]->ReadFile("xtreeTim.dat",xstr);

        xtree[f]->SetBranchAddress("evnum",&evnum[f]);
        xtree[f]->SetBranchAddress("Q2",&Q2[f]);
        xtree[f]->SetBranchAddress("x",&x[f]);
        xtree[f]->SetBranchAddress("W",&W[f]);
        xtree[f]->SetBranchAddress("y",&y[f]);
        xtree[f]->SetBranchAddress("eleP",&eleP[f]);
        xtree[f]->SetBranchAddress("Mh",&Mh[f]);
        xtree[f]->SetBranchAddress("xF",&xF[f]);
        xtree[f]->SetBranchAddress("PhPerp",&PhPerp[f]);
        xtree[f]->SetBranchAddress("theta",&theta[f]);
        xtree[f]->SetBranchAddress("PhiH",&PhiH[f]);
        xtree[f]->SetBranchAddress("PhiR",&PhiR[f]);
        for(h=0; h<nHad; h++) {
          xtree[f]->SetBranchAddress( TString(hadN[h]+"P"), &hadP[f][h] );
        };

        break;

      case kHarut:
        compareTitle[f] = "Harut -- dihad1000.dis.0000.nrad.dat.evio.hipo2.txt";
        compareName[f] = "harut";
        gROOT->ProcessLine(".! python formatHarutFile.py");
        xstr = "evnum/I";
        for(h=0; h<nHad; h++) {
          xstr += ":"+hadN[h]+"E/F";
          xstr += ":"+hadN[h]+"Theta/F";
          xstr += ":"+hadN[h]+"Phi/F";
          xstr += ":"+hadN[h]+"Ptq/F";
          xstr += ":"+hadN[h]+"XF/F";
        };
        xstr += ":PhPerp/F:PhiH/F:PhiR1/F:PhiR2/F";
        printf("xtree[%d] branches: %s\n",f,xstr.Data());
        xtree[f]->ReadFile("xtreeHarut.dat",xstr);

        xtree[f]->SetBranchAddress("evnum",&evnum[f]);
        xtree[f]->SetBranchAddress("PhPerp",&PhPerp[f]);
        xtree[f]->SetBranchAddress("PhiH",&PhiH[f]);
        xtree[f]->SetBranchAddress("PhiR2",&PhiR[f]);
        for(h=0; h<nHad; h++) {
          xtree[f]->SetBranchAddress( TString(hadN[h]+"E"), &hadE[f][h] );
          xtree[f]->SetBranchAddress( TString(hadN[h]+"Ptq"), &hadPtq[f][h] );
          xtree[f]->SetBranchAddress( TString(hadN[h]+"Theta"), &hadTheta[f][h] );
          xtree[f]->SetBranchAddress( TString(hadN[h]+"Phi"), &hadPhi[f][h] );
          xtree[f]->SetBranchAddress( TString(hadN[h]+"XF"), &hadXF[f][h] );
        };

        break;

      case kOrlando:
        // Orlando's data
        compareTitle[f] = "Orlando -- dataOS_HA_mc.txt";
        compareName[f] = "orlando";
        ReadOrlandoFormat("xfiles/dataOS_HA_mc.txt",f);
        break;

      case kHarutOS:
        // Harut's data, reformatted by Orlando
        compareTitle[f] = "Harut (formatted by Orlando) -- dataHA.txt";
        compareName[f] = "harutOS";
        ReadOrlandoFormat("xfiles/dataHA.txt",f);
        break;


      case kSimpleC:
        compareTitle[f] = "Chris (simple tree from Clas12Tool) -- simpleTree.root";
        compareName[f] = "chrisSimpleC";

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
        compareTitle[f] = "Chris (simple tree from java tools) -- javaOut.root";
        compareName[f] = "chrisSimpleJava";
        
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
      
      // using EventTree these cases; no need to set branch addresses
      case kAnalysis: 
        compareTitle[f] = "Chris (REC::Particle bank) -- outrootCC/rec.root";
        compareName[f] = "chrisRec";
        break;
      case kAnalysisLund:
        compareTitle[f] = "Chris (MC::Lund bank) -- outrootCC/lund.root";
        compareName[f] = "chrisLund";
        break;

      default: 
        fprintf(stderr,"ERROR: unknown xcheck[%d]=%d\n",f,xcheck[f]);
        exit(0);

    };
  };

  for(f=0; f<NF; f++) {
    if(xcheck[f]==kAnalysis || xcheck[f]==kAnalysisLund) ENT[f] = ev[f]->ENT;
    else ENT[f] = xtree[f]->GetEntries();
  };




  // build hash table for xtree[1]
  // -- later we will loop through xtree[0], searching for each event's hash value in
  //    this hash table in order to find a matching event in xtree[1]
  // -- this hash function now just returns the event number; before when our event
  //    numbers were not matching, I was instead using a hash function to search for
  //    matching events
  Int_t hashVal;
  std::map<Int_t,Int_t> hashMap; // event hash -> xtree[1] index
  std::map<Int_t,Int_t>::iterator hashIter;

  for(int xi=0; xi<ENT[1]; xi++) {

    GetTreeVars(1,xi);

    // hash function
    switch(xcheck[1]) {
      default: hashVal = evnum[1];
    };

    hashMap.insert(std::pair<Int_t,Int_t>(hashVal,xi));
  };



  // define output file
  TString outdat = "compare_" + compareName[0] + "_" + compareName[1] + ".dat";
  gSystem->RedirectOutput(outdat,"w");
  for(f=0; f<NF; f++) printf("xtree%d:  %s\n",f,compareTitle[f].Data());
  printf("\n");
  gSystem->RedirectOutput(0);

  Bool_t extraCut,evCut;
  

  // loop through xtree[0]
  for(int i=0; i<ENT[0]; i++) {

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
      PhPerp[f] = -10000;
      eleP[f] = -10000;
      for(h=0; h<nHad; h++) {
        hadE[f][h] = -10000;
        hadP[f][h] = -10000;
        hadPt[f][h] = -10000;
        hadPtq[f][h] = -10000;
        hadTheta[f][h] = -10000;
        hadPhi[f][h] = -10000;
        hadPx[f][h] = -10000;
        hadPy[f][h] = -10000;
        hadPz[f][h] = -10000;
        hadZ[f][h] = -10000;
        hadXF[f][h] = -10000;
      };
    };

    // fill xtree[0] variables
    GetTreeVars(0,i);

    //ev->GetEvent(i); // evtr loop


    evCut = true; // TODO
    //if(xcheck[0]==kAnalysis) evCut = ev->cutCrossCheck;
    //else... 
    if(evCut) {

      // hash xtree[0]'s event
      switch(xcheck[0]) {
        default: hashVal = evnum[0];
      };


      // check if xtree[1] has a matching hash value
      hashIter = hashMap.find(hashVal);
      if(hashIter!=hashMap.end()) {

        // set xtree[1] variables to the matching event's
        GetTreeVars(1,hashIter->second);

        // extra requirement to improve event matching
        if( xcheck[0]==kHarutOS || xcheck[1]==kHarutOS || 
            xcheck[0]==kHarut || xcheck[1]==kHarut ||
            xcheck[0]==kOrlando || xcheck[1]==kOrlando
        ) {
          extraCut = fabs(hadE[0][iP]-hadE[1][iP]) < 0.1 &&
                     fabs(hadE[0][iM]-hadE[1][iM]) < 0.1; //////////////////!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        }
        else extraCut = true;


        if(extraCut) {

          // print comparisons
          // -- if a variable is not set (i.e., set to -10000), its comparison will not
          //    be printed
          gSystem->RedirectOutput(outdat,"a");

          printf("EVENT#");
          for(f=0; f<NF; f++) printf("  xtree%d: %d",f,evnum[f]); printf("\n");
          //printf("HASH  xtree0: %d  xtree1: %d\n",hashVal,hashIter->first);
          printf("%12s %12s %12s %12s %12s\n",
            "evnum","var",compareName[0].Data(),compareName[1].Data(),"diff");
          printf("%12s %12s %12s %12s %12s\n",sep,sep,sep,sep,sep);


          for(h=0; h<nHad; h++) {
            PrintCompare( TString(hadN[h]+"E"), hadE[0][h], hadE[1][h] );
            PrintCompare( TString(hadN[h]+"P"), hadP[0][h], hadP[1][h] );
            PrintCompare( TString(hadN[h]+"Pt"), hadPt[0][h], hadPt[1][h] );
            PrintCompare( TString(hadN[h]+"Ptq"), hadPtq[0][h], hadPtq[1][h] );
            PrintCompare( TString(hadN[h]+"Px"), hadPx[0][h], hadPx[1][h] );
            PrintCompare( TString(hadN[h]+"Py"), hadPy[0][h], hadPy[1][h] );
            PrintCompare( TString(hadN[h]+"Pz"), hadPz[0][h], hadPz[1][h] );
            PrintCompare( TString(hadN[h]+"Theta"), hadTheta[0][h], hadTheta[1][h] );
            PrintCompare( TString(hadN[h]+"Phi"), hadPhi[0][h], hadPhi[1][h] );
            PrintCompare( TString(hadN[h]+"Z"), hadZ[0][h], hadZ[1][h] );
            PrintCompare( TString(hadN[h]+"XF"), hadXF[0][h], hadXF[1][h] );
          };

          PrintCompare( "eleP", eleP[0], eleP[1] );
          PrintCompare( "Q2", Q2[0], Q2[1] );
          PrintCompare( "W", W[0], W[1] );
          PrintCompare( "x", x[0], x[1] );
          PrintCompare( "y", y[0], y[1] );
          PrintCompare( "Mh", Mh[0], Mh[1] );
          PrintCompare( "xF", xF[0], xF[1] );
          PrintCompare( "PhPerp", PhPerp[0], PhPerp[1] );
          PrintCompare( "theta", theta[0], theta[1] );
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

  if(name.Contains("Phi")) {
    // if it's an angle, ensure it's in proper range
    val0 = Tools::AdjAngle(val0);
    val1 = Tools::AdjAngle(val1);
    diff = Tools::AdjAngle( val0 - val1 );
  } else {
    // otherwise just compare the values
    diff = val0 - val1;
  };

  TString suffix = fabs(diff) > 1e-4 ? "  <- disagreement" : "";

  // print comparison
  printf("%12d %12s %12.5f %12.5f %12.5f%s\n",
    evnum[0],name.Data(),val0,val1,diff,suffix.Data());

};


// set local variables to xtree (or EventTree) variables
void GetTreeVars(Int_t ff, Int_t ii) {
  if(xcheck[ff]==kAnalysis || xcheck[ff]==kAnalysisLund) GetEventTreeVars(ff,ii);
  else xtree[ff]->GetEntry(ii);
};

// set local variables to EventTree variables
void GetEventTreeVars(Int_t ff, Int_t ii) {
  ev[ff]->GetEvent(ii);
  if(ev[ff]->pairType == EncodePairType(kPip,kPim)) {
    evnum[ff] = ev[ff]->evnum;
    Q2[ff] = ev[ff]->Q2;
    W[ff] = ev[ff]->W;
    x[ff] = ev[ff]->x;
    y[ff] = ev[ff]->y;
    Mh[ff] = ev[ff]->Mh;
    xF[ff] = ev[ff]->xF;
    PhPerp[ff] = ev[ff]->PhPerp;
    theta[ff] = ev[ff]->theta;
    PhiH[ff] = ev[ff]->PhiH;
    PhiR[ff] = ev[ff]->PhiR;
    eleP[ff] = ev[ff]->eleP;
    for(int hh=0; hh<nHad; hh++) {
      hadE[ff][hh] = ev[ff]->hadE[hh];
      hadP[ff][hh] = ev[ff]->hadP[hh];
      hadPt[ff][hh] = ev[ff]->hadPt[hh];
      hadPtq[ff][hh] = ev[ff]->hadPtq[hh];
      hadTheta[ff][hh] = Tools::EtaToTheta(ev[ff]->hadEta[hh]);
      hadPhi[ff][hh] = ev[ff]->hadPhi[hh];
      hadZ[ff][hh] = ev[ff]->Z[hh];
      hadXF[ff][hh] = ev[ff]->hadXF[hh];
    };
  };
};


// read Orlando's format table 
// (Orlando formatted Harut's data into a table, which has the same
// format as Orlando's table, so here's a method that can read either one
// of them, so that we don't have 2 copies of the same code
void ReadOrlandoFormat(TString datname, Int_t ff) {

  // strip header row (easier/lazier than editting header and using TTree::ReadFile)
  TString datoutname = datname;
  datoutname.ReplaceAll("xfiles/","xtree");
  datoutname.ReplaceAll(".txt",".dat");
  gROOT->ProcessLine(TString(".! tail -n +2 " + datname + " > " + datoutname));

  // set branch names
  xstr = "evnum/I";
  xstr += ":pipE/F:pipTheta/F:pipPhi/F:pipPhiH/F:pipPtq/F";
  xstr += ":pipZ/F:pipMx/F:pipXF/F:pipEtaCM/F:pipEtaBreit/F";
  xstr += ":pimE/F:pimTheta/F:pimPhi/F:pimPhiH/F:pimPtq/F";
  xstr += ":pimZ/F:pimMx/F:pimXF/F:pimEtaCM/F:pimEtaBreit/F";
  xstr += ":PhiH/F:PhPerp/F:PhiR1/F:PhiR2/F";

  printf("xtree[%d] branches: %s\n",ff,xstr.Data());
  xtree[ff]->ReadFile(datoutname,xstr);

  xtree[ff]->SetBranchAddress("evnum",&evnum[ff]);
  for(int hh=0; hh<nHad; hh++) {
    xtree[ff]->SetBranchAddress( TString(hadN[hh]+"E"), &hadE[ff][hh] );
    xtree[ff]->SetBranchAddress( TString(hadN[hh]+"Theta"), &hadTheta[ff][hh] );
    xtree[ff]->SetBranchAddress( TString(hadN[hh]+"Phi"), &hadPhi[ff][hh] );
    xtree[ff]->SetBranchAddress( TString(hadN[hh]+"Ptq"), &hadPtq[ff][hh] );
    xtree[ff]->SetBranchAddress( TString(hadN[hh]+"Z"), &hadZ[ff][hh] );
    xtree[ff]->SetBranchAddress( TString(hadN[hh]+"XF"), &hadXF[ff][hh] );
  };
  xtree[ff]->SetBranchAddress("PhPerp",&PhPerp[ff]);
  xtree[ff]->SetBranchAddress("PhiH",&PhiH[ff]);
  //xtree[ff]->SetBranchAddress("PhiR1",&PhiR[ff]);
  xtree[ff]->SetBranchAddress("PhiR2",&PhiR[ff]); // preferred?
};
