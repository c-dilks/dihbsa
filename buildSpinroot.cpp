// builds spinroot files, which contain the necessary data structures
// used by the asymmetry fitting code

#include <cstdlib>
#include <iostream>
#include <vector>
#include <map>

// ROOT
#include "TFile.h"
#include "TString.h"
#include "TRegexp.h"

// DihBsa
#include "Constants.h"
#include "EventTree.h"
#include "Binning.h"
#include "Asymmetry.h"


// argument variables
TString inputData;
Int_t pairType;
Int_t whichModulation;
Int_t ivType;
Int_t whichHelicityMC;

// subroutines
void SetDefaultArgs();
int PrintUsage();

// global variables
Binning * BS;
Asymmetry * A;
EventTree * ev;


//////////////////////////////////////


int main(int argc, char** argv) {

  //gDebug = 2; // use to debug streaming problems
  
  // read options
  SetDefaultArgs();
  int opt;
  enum inputType_enum {iFile,iDir};
  Int_t inputType = -1;
  while( (opt=getopt(argc,argv,"f:d:p:m:i:h:")) != -1 ) {
    switch(opt) {
      case 'f': /* input file */
        if(inputType>=0) return PrintUsage();
        inputData = optarg;
        inputType = iFile;
        break;
      case 'd': /* input directory */
        if(inputType>=0) return PrintUsage();
        inputData = optarg;
        inputType = iDir;
        break;
      case 'p': /* pair type (hexadecimal number) */
        pairType = (Int_t) strtof(optarg,NULL);
        break;
      case 'm': /* azimuthal modulation for asymmetry */
        whichModulation = (Int_t) strtof(optarg,NULL);
        break;
      case 'i': /* independent variables */
        ivType = (Int_t) strtof(optarg,NULL);
        break;
      case 'h': /* which helicityMC */
        whichHelicityMC = (Int_t) strtof(optarg,NULL);
        break;
      default: return PrintUsage();
    };
  };
  if(inputType!=iFile && inputType!=iDir) {
    fprintf(stderr,"ERROR: must specify input file or directory\n");
    return PrintUsage();
  };
  if(pairType==0x99) {
    fprintf(stderr,"ERROR: pi0 pi0 channel is under construction\n");
    return 0;
  };

  // print arguments' values
  printf("inputData = %s\n",inputData.Data());
  printf("pairType = 0x%x\n",pairType);
  printf("whichModulation = %d\n",whichModulation);
  printf("ivType = %d\n",ivType);
  printf("whichHelicityMC = %d\n",whichHelicityMC);
  printf("\n");


  // set binning scheme
  BS = new Binning(pairType);
  BS->AsymModulation = whichModulation;
  Bool_t schemeSuccess = BS->SetScheme(ivType);
  if(!schemeSuccess) {
    fprintf(stderr,"ERROR: Binning::SetScheme failed\n");
    return 0;
  };


  // instantiate EventTree 
  // (use 1 file if inputType==iFile, or all root files in inputData if inputType==iDir)
  ev = new EventTree(inputData+(inputType==iDir?"/*.root":""),pairType);


  // get modulation name for 1-amp fit
  A = new Asymmetry(BS);
  TString modN = A->ModulationName;
  printf("--> 1-amp fit will be for %s modulation\n\n",modN.Data());


  // instantiate spinroot file
  TString spinrootFileN;
  if(inputType==iDir) {
    spinrootFileN = "spinroot/spin";
    spinrootFileN += "__" + PairName(pairType) + "_";
    for(int d=0; d<BS->dimensions; d++) spinrootFileN += "_" + BS->GetIVname(d);
    spinrootFileN += "__" + modN;
    spinrootFileN += ".root";
  } else if(inputType==iFile) {
    spinrootFileN = inputData;
    spinrootFileN(TRegexp("^.*/")) = "spinroot/spin.";
  };
  printf("\nCREATING OUTPUT FILE = %s\n\n",spinrootFileN.Data());
  TFile * spinrootFile = new TFile(spinrootFileN,"RECREATE");


  // instantiate Asymmetry objects, and
  // build map of 3-digit bin number -> Asymmetry object
  std::map<Int_t, Asymmetry*> asymMap;
  for(Int_t bn : BS->binVec) {
    A = new Asymmetry(BS,bn);
    if(A->success) asymMap.insert(std::pair<Int_t, Asymmetry*>(bn,A));
    else return 0;
  };
        

  //-----------------------------------------------------
  // EVENT LOOP  
  //-----------------------------------------------------
  Bool_t eventAdded;
  ev->whichHelicityMC = whichHelicityMC;

  printf("begin loop through %lld events...\n",ev->ENT);
  for(int i=0; i<ev->ENT; i++) {

    ev->GetEvent(i);
    if(ev->Valid()) {

      // add event to Asymmetry, where Asymmetry::AddEvent() checks the bin,
      // and fills plots if it's the correct bin; if not, AddEvent() does nothing
      for(Int_t bn : BS->binVec) {
        A = asymMap.at(bn);
        eventAdded = A->AddEvent(ev);
        //if(eventAdded && A->debug) ev->PrintEvent();
      };

    };
  }; // eo EVENT LOOP


  // write out to spinroot file
  spinrootFile->cd();
  BS->Write("BS");
  for(Int_t bn : BS->binVec) {
    A = asymMap.at(bn);
    A->StreamData(spinrootFile);
  };


  // stream objects directly to spinroot file
  // - this doesn't work, unfortunately. The Asymmetry objects get written, and have a
  //   data size that indicates histograms/data structures are being written, however
  //   upon trying to access anything in the TFile, it immediately seg faults
  // - alternative implementation above writes out only pertinent members of Asymmetry
  //   instead of the full Asymmetry object
  /*
  TString Aname;
  spinrootFile->cd();
  BS->Write("BS");
  for(Int_t bn : BS->binVec) {
    A = asymMap.at(bn);
    Aname = "A" + A->binN;
    printf("write %s\n",Aname.Data());
    A->PrintSettings();
    //A->rfData->Write(Aname);
    A->Write(Aname);
  };
  */


  // close spinroot file
  spinrootFile->Close();
  printf("\n%s written\n\n",spinrootFileN.Data());
};

//////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////


// set default arguments
void SetDefaultArgs() {
  inputData = "";
  pairType = EncodePairType(kPip,kPim);
  whichModulation = Asymmetry::weightSinPhiHR;
  ivType = Binning::vM + 1;
  whichHelicityMC = 0;
};


// help printout
int PrintUsage() {

  SetDefaultArgs();
  BS = new Binning(EncodePairType(kPip,kPim));
  fprintf(stderr,"\nUSAGE: buildSpinroot.exe [-f or -d input_data ] [options...]\n\n");

  printf("INPUT DATA:\n");
  printf(" -f\tsingle ROOT file\n");
  printf(" -d\tdirectory of ROOT files\n");
  printf(" NOTE: specify input with either -f or -d, but not both\n");
  printf("\n");

  printf("OPTIONS:\n");

  printf(" -p\tpair type, specified as a hexadecimal number\n");
  printf("   \trun PrintEnumerators.C for notation\n");
  printf("   \tdefault = 0x%x (%s)\n\n",pairType,PairTitle(pairType).Data());

  printf(" -m\tazimuthal modulation for asymmetry linear fit\n");
  BS = new Binning(pairType);
  for(int m=0; m<Asymmetry::nMod; m++) {
    BS->AsymModulation = m;
    A = new Asymmetry(BS);
    printf("   \t %d = %s =  %s\n",m,
        (A->ModulationName).Data(),
        (A->ModulationTitle).Data()
        );
  };
  printf("   \tdefault = %d\n\n",whichModulation);

  printf(" -i\tindependent variable specifier: 1, 2, or 3-digit number which\n");
  printf("   \tspecifies the independent variables that asymmetries will be\n");
  printf("   \tplotted against. The number of digits will be the number of\n");
  printf("   \tdimensions in the multi-dimensional binning\n");
  printf("   \t* the allowed digits are:\n");
  for(int i=0; i<Binning::nIV; i++) {
    printf("   \t  %d = %s\n",i+1,(BS->IVtitle[i]).Data());
  };
  printf("   \tdefault = %d\n\n",ivType);

  printf(" -g\t (for MC) - select which helicityMC to use\n\n");

  return 0;
};

