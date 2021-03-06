#include <cstdlib>
#include <iostream>
#include <vector>
#include <map>

// ROOT
#include "TFile.h"
#include "TTree.h"
#include "TChain.h"
#include "TString.h"
#include "TMath.h"
#include "TSystem.h"
#include "TRegexp.h"
#include "TH1.h"
#include "TH2.h"
#include "TF1.h"
#include "TCanvas.h"
#include "TLine.h"
#include "TStyle.h"
#include "TMultiGraph.h"
#include "TSystemDirectory.h"
#include "TSystemFile.h"
#include "TList.h"
#include "TCollection.h" // for TIter

// DihBsa
#include "Constants.h"
#include "DIS.h"
#include "Trajectory.h"
#include "Dihadron.h"
#include "EventTree.h"
#include "Binning.h"
#include "Asymmetry.h"

// subroutines
Int_t HashBinNum(Int_t bin0, Int_t bin1=-1, Int_t bin2=-1);
void DrawKinDepGraph(TGraphErrors * g_, Binning * B_, Int_t v_);
void DrawSimpleGraph(TGraphErrors * g_, Binning * B_, Int_t v, Bool_t setRange_=true);
void DrawAsymGr(TGraphErrors * g_);
void DrawAsymGr2(TGraph2DErrors * g_);
void SetCloneName(TH1 * clone_);
TGraphErrors * ShiftGraph(TGraphErrors * gr, Int_t nShift);
int PrintUsage();
void SetDefaultArgs();

// argument variables
TString inputData;
Int_t pairType;
Int_t whichModulation;
Int_t ivType;
Int_t flowControl;
Int_t whichHelicityMC;

// other global variables
Int_t dimensions;
Int_t N_AMP,N_D;
Int_t whichHad[2];
TString dihTitle, dihName;
Binning * BS;
Asymmetry * A;

enum flowEnum { 
  fSerial,
  fSerialRenamed,
  fParallelFill,
  fParallelCat,
  fParallelCalc
}; // for flowControl setting


//////////////////////////////////////


int main(int argc, char** argv) {

  gStyle->SetOptFit(1);
  //gDebug = 2; // use to debug streaming problems

  SetDefaultArgs();


  // read options
  int opt;
  Int_t inputType = 0; // 1=singleFile, 2=directory
  while( (opt=getopt(argc,argv,"f:d:p:m:i:c:h:")) != -1 ) {
    switch(opt) {
      case 'f': /* input file */
        if(inputType>0) return PrintUsage();
        inputData = optarg;
        inputType = 1;
        break;
      case 'd': /* input directory */
        if(inputType>0) return PrintUsage();
        inputData = optarg;
        inputType = 2;
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
      case 'c': /* flow control */
        flowControl = (Int_t) strtof(optarg,NULL);
        break;
      case 'h': /* which helicityMC */
        whichHelicityMC = (Int_t) strtof(optarg,NULL);
        break;
      default: return PrintUsage();
    };
  };

  if(inputType==0 && flowControl!=fParallelCat && flowControl!=fParallelCalc) {
    fprintf(stderr,"ERROR: must specify input file or directory\n");
    return PrintUsage();
  };



  // print arguments' values
  printf("inputData = %s\n",inputData.Data());
  printf("pairType = 0x%x\n",pairType);
  printf("whichModulation = %d\n",whichModulation);
  printf("ivType = %d\n",ivType);
  printf("flowControl = %d\n",flowControl);
  printf("whichHelicityMC = %d\n",whichHelicityMC);
  printf("\n");

  // pi0pi0 not yet functional TODO
  if(pairType==0x99) {
    fprintf(stderr,"ERROR: pi0 pi0 channel is under construction\n");
    return 0;
  };


  // set dihadron name / title
  DecodePairType(pairType,whichHad[qA],whichHad[qB]);
  dihTitle = PairTitle(whichHad[qA],whichHad[qB]);
  dihName = PairName(whichHad[qA],whichHad[qB]);


  // determine number of dimensions (for multi-dimensional binning)
  dimensions=1;
  if(ivType>=10) dimensions=2;
  if(ivType>=100) dimensions=3;
  if(ivType>=1000) {
    fprintf(stderr,"ERROR: ivType has too many digits\n");
    return 0;
  };



  // ivType:
  // -- if dimensions==1, all asymmetries will be plotted against one independent
  //    variable (IV); in this case, ivType is that IV, according to enumerators in
  //    Binning (ivEnum)
  // -- if dimensions==2, asymmetries are plotted for two IVs. For this one, ivType is
  //    understood as a 2-digit number: the first digit is IV0, and the second is IV1.
  //    The asymmetries will be plotted vs. IV0, for bins in IV1
  // -- if dimensions==3, we have 3 IVs and ivType is a 3-digit number, representing
  //    IV0, IV1, and IV2. Aymmetries are plotted vs IV0 for bins in (IV1,IV2) pairs

  // read ivType and convert it to IV enumerators
  Int_t ivVar[3] = {-1,-1,-1};
  switch(dimensions) {
    case 1:
      ivVar[0] = ivType - 1;
      break;
    case 2:
      ivVar[0] = ivType / 10 - 1; 
      ivVar[1] = ivType % 10 - 1;
      break;
    case 3:
      ivVar[0] = ivType / 100 - 1;
      ivVar[1] = ( ivType / 10 ) % 10 - 1;
      ivVar[2] = ivType % 10 - 1;
      break;
    default:
      fprintf(stderr,"ERROR: bad number of dimensions\n");
      return 0;
  };


  // check IV enumerators and get number of bins for each IV
  BS = new Binning(pairType); // instantiate binning scheme 

  // DNP2019 shorctut: interested in z-dependence above and below rho mass
  // - override binning scheme, if we requested z-dependence in mass bins
  Int_t vvm = Binning::vM;
  if(ivType==32) {
    BS->nBins[vvm] = -1;
    BS->bound[vvm].clear();
    BS->AddBinBound(vvm,BS->minIV[vvm]);
    BS->AddBinBound(vvm,0.77);
    BS->AddBinBound(vvm,BS->maxIV[vvm]);
    BS->PrintBinBounds();
  };
  //// END shortcut

  Int_t NB[3]; // # of bins for each IV
  for(int d=0; d<dimensions; d++) {
    if(!(BS->ValidIV(ivVar[d]))) {
      printf("ivVar[%d] = %d\n",d,ivVar[d]);
      fprintf(stderr,"ERROR: this IV is unknown\n");
      return 0;
    };
    NB[d] = BS->nBins[ivVar[d]];
  };


  // print which IV will be analyzed
  A = new Asymmetry(BS,whichModulation,UNDEF);
  TString modN = A->ModulationName;
  printf("--------> Analysing %s asymmetries vs. %s ",
      modN.Data(),(BS->IVname[ivVar[0]]).Data());
  if(dimensions>=2)
    printf("in bins of %s ",(BS->IVname[ivVar[1]]).Data());
  if(dimensions>=3)
    printf("and %s ",(BS->IVname[ivVar[2]]).Data());
  printf("\n\n");


  // instantiate EventTree
  EventTree * ev;
  if(flowControl==fSerial || flowControl==fSerialRenamed) {
    if(inputType==2) 
      ev = new EventTree(TString(inputData+"/*.root"),pairType);
    else {
      fprintf(stderr,"ERROR: flowControl is serial, inputData must be a directory\n");
      return 0;
    };
  } else if(flowControl==fParallelFill) {
    if(inputType==1) 
      ev = new EventTree(inputData,pairType);
    else {
      fprintf(stderr,"ERROR: flowControl is parallel, inputData must be a file\n");
      return 0;
    };
  };


  // set output file names 
  TString outName;
  TFile *resultFile, *spinrootFile;
  if(flowControl==fSerial || flowControl==fSerialRenamed) {
    outName = "spin";
    if(flowControl==fSerialRenamed) {
      outName = "spinout/" + outName;
      outName += "__" + dihName + "_";
      for(int d=0; d<dimensions; d++) outName += "_" + BS->IVname[ivVar[d]];
      outName += "__" + modN;
    };
    outName += ".root";
    resultFile = new TFile(outName,"RECREATE");
  }
  else if(flowControl==fParallelFill) {
    outName = inputData;
    outName(TRegexp("^.*/")) = "spinroot/spin.";
    printf("outName = %s\n",outName.Data());
    spinrootFile = new TFile(outName,"RECREATE");
  }
  else if(flowControl==fParallelCat) {
    spinrootFile = new TFile("spinrootCat.root","RECREATE");
  }
  else if(flowControl==fParallelCalc) {
    resultFile = new TFile("spinFinal.root","RECREATE");
  }
  else {
    fprintf(stderr,"ERROR: bad flowControl (%d)\n",flowControl);
    return 0;
  };
  printf("\nCREATING OUTPUT FILE = %s\n\n",outName.Data());




  // MAPS:
  // - each IV bin is assigned a 3-digit bin number
  // - the maps map that 3-digit bin number to an object (Asymmetry object, plot, etc)
  // - use HashBinNum(...) to convert a set of bin numbers to the 3-digit bin number
  std::vector<Int_t> binVec; // vector of 3-digit bin numbers, used for looping
  // over all possible bins
  std::map<Int_t, Int_t> binVecMap[3]; // map 3-digit bin number back to each
  // single-digit bin number [dimension]
  std::map<Int_t, Asymmetry*> asymMap; // map to Asymmetry pointers
  std::map<Int_t, TGraphErrors*> kindepMap; // map to kindep plots
  std::map<Int_t, TGraphErrors*> RFkindepMap[Asymmetry::nAmp]; 
  std::map<Int_t, TMultiGraph*> multiMap; // map to multiGr
  std::map<Int_t, TGraphErrors*> chindfMap; // map to chindfGr
  std::map<Int_t, TGraphErrors*> rellumMap; // map to rellumGr
  // map to MLM amplitudes [for each param]

  TGraphErrors * kindepGr;
  TGraphErrors * RFkindepGr[Asymmetry::nAmp];
  TGraphErrors * chindfGr;
  TGraphErrors * rellumGr;
  TMultiGraph * multiGr;

  TString grTitle,grTitleSuffix,grName,grNameSuffix;
  Int_t binn[3];
  Bool_t first;



  // build binVec and instantiate Asymmetry objects
  Int_t binNum;
  if(dimensions == 1) {
    for(int b=0; b<NB[0]; b++) {

      binNum = HashBinNum(b);
      binVec.push_back(binNum);
      binVecMap[0].insert(std::pair<Int_t,Int_t>(binNum,b));
      binVecMap[1].insert(std::pair<Int_t,Int_t>(binNum,-1));
      binVecMap[2].insert(std::pair<Int_t,Int_t>(binNum,-1));

      A = new Asymmetry(BS, whichModulation, 1, ivVar[0], b);
      if(A->success) asymMap.insert(std::pair<Int_t, Asymmetry*>(binNum,A));
      else return 0;

    };
  }
  else if(dimensions == 2) {
    for(int b1=0; b1<NB[1]; b1++) {
      for(int b0=0; b0<NB[0]; b0++) {

        binNum = HashBinNum(b0,b1);
        binVec.push_back(binNum);
        binVecMap[0].insert(std::pair<Int_t,Int_t>(binNum,b0));
        binVecMap[1].insert(std::pair<Int_t,Int_t>(binNum,b1));
        binVecMap[2].insert(std::pair<Int_t,Int_t>(binNum,-1));

        A = new Asymmetry(BS, whichModulation, 2, 
            ivVar[0], b0,
            ivVar[1], b1
            );
        if(A->success) asymMap.insert(std::pair<Int_t, Asymmetry*>(binNum,A));
        else return 0;

      };
    };
  }
  else if(dimensions == 3) {
    for(int b2=0; b2<NB[2]; b2++) {
      for(int b1=0; b1<NB[1]; b1++) {
        for(int b0=0; b0<NB[0]; b0++) {

          binNum = HashBinNum(b0,b1,b2);
          binVec.push_back(binNum);
          binVecMap[0].insert(std::pair<Int_t,Int_t>(binNum,b0));
          binVecMap[1].insert(std::pair<Int_t,Int_t>(binNum,b1));
          binVecMap[2].insert(std::pair<Int_t,Int_t>(binNum,b2));

          A = new Asymmetry(BS, whichModulation, 3, 
              ivVar[0], b0,
              ivVar[1], b1,
              ivVar[2], b2
              );
          if(A->success) asymMap.insert(std::pair<Int_t, Asymmetry*>(binNum,A));
          else return 0;

        };
      };
    };
  };


  // instantiate all graphs etc. and insert them in the map
  first = true;
  for(Int_t bn : binVec) {

    for(int d=0; d<3; d++) binn[d]=binVecMap[d].at(bn);
     
    A = asymMap.at(bn);

    // get bin numbers from 3-digit bin number; get number of RooFit params
    if(first) { 
      N_AMP = A->nAmpUsed;
      N_D = A->nDparamUsed;
    };
    first = false;


    // set title and name suffixes
    switch(dimensions) {
      case 1:
        grTitleSuffix =
          BS->IVtitle[ivVar[0]] + ";" + BS->IVtitle[ivVar[0]];
        grNameSuffix = BS->IVname[ivVar[0]];
        break;
      case 2:
        grTitleSuffix = BS->IVtitle[ivVar[0]] + " :: " +
          BS->GetBoundStr(ivVar[1],binn[1]) + ";" +
          BS->IVtitle[ivVar[0]];
        grNameSuffix = Form("%s_bin_%s%d",
            (BS->IVname[ivVar[0]]).Data(),
            (BS->IVname[ivVar[1]]).Data(), binn[1]
            );
        break;
      case 3:
        grTitleSuffix = BS->IVtitle[ivVar[0]] + " :: " +
          BS->GetBoundStr(ivVar[1],binn[1]) + ", " +
          BS->GetBoundStr(ivVar[2],binn[2]) + ";" +
          BS->IVtitle[ivVar[0]];
        grNameSuffix = Form("%s_bin_%s%d_%s%d",
            (BS->IVname[ivVar[0]]).Data(),
            (BS->IVname[ivVar[1]]).Data(), binn[1],
            (BS->IVname[ivVar[2]]).Data(), binn[2]
            );
        break;
    };


    // instantiate graphs; only needs to be done for each IV1 and IV2 bin, since the
    // horizontal axis of these graphs are all IV0 (hence the requirement binn[0]==0)
    if(binn[0]==0) {

      // instantiate kindep graph, for linear fit ("l.f.") result
      grTitle = dihTitle + " A_{LU}[" + A->ModulationTitle + "]_{l.f.} " + 
        " vs. " + grTitleSuffix;
      grName = "kindep_" + grNameSuffix;
      kindepGr = new TGraphErrors();
      kindepGr->SetName(grName);
      kindepGr->SetTitle(grTitle);

      // instantiate multiGraph, for plotting kindep graphs together
      grTitle = dihTitle + " A_{LU}[" + A->ModulationTitle + "] " + 
        " vs. " + grTitleSuffix;
      grName = "multiGr_" + grNameSuffix;
      multiGr = new TMultiGraph();
      multiGr->SetName(grName);
      multiGr->SetTitle(grTitle);

      // instantiate kindep graphs for maximum likelihood method (m.l.m.) result
      // for each fit parameter
      for(int aa=0; aa<N_AMP; aa++) {
        grTitle = dihTitle + " " + TString(A->rfA[aa]->GetTitle()) + "_{m.l.m.} " +
          " vs. " + grTitleSuffix;
        grName = "RF_A" + TString::Itoa(aa,10) + "_kindep_" + grNameSuffix;
        RFkindepGr[aa] = new TGraphErrors();
        RFkindepGr[aa]->SetName(grName);
        RFkindepGr[aa]->SetTitle(grTitle);
      };

      // instantiate chi2/ndf graphs
      grTitle = "#chi^{2}/NDF of " +
        dihTitle + " A_{LU}[" + A->ModulationTitle + "]_{l.f.} " + 
        " vs. " + grTitleSuffix;
      grName = "chindf_" + grNameSuffix;
      chindfGr = new TGraphErrors();
      chindfGr->SetName(grName);
      chindfGr->SetTitle(grTitle);

      // instantiate relative luminosity graphs
      grTitle = "relative luminosity vs. " + grTitleSuffix;
      grName = "rellum_" + grNameSuffix;
      rellumGr = new TGraphErrors();
      rellumGr->SetName(grName);
      rellumGr->SetTitle(grTitle);

    }; // eo if(binn[0]==0)


    // insert objects into maps (note: these are many-to-one, i.e., several
    // bin numbers will map to the same pointer)
    kindepMap.insert(std::pair<Int_t,TGraphErrors*>(bn,kindepGr));
    multiMap.insert(std::pair<Int_t,TMultiGraph*>(bn,multiGr));
    chindfMap.insert(std::pair<Int_t,TGraphErrors*>(bn,chindfGr));
    rellumMap.insert(std::pair<Int_t,TGraphErrors*>(bn,rellumGr));
    for(int aa=0; aa<N_AMP; aa++) {
      RFkindepMap[aa].insert(std::pair<Int_t,TGraphErrors*>(bn,RFkindepGr[aa]) );
    };



  }; // eo binVec loop



  // overall summary plots
  TH1F * chisqDist = new TH1F("chisqDist",
      "#chi^{2} distribution (from linear fit results)",100,0,20);



  //-----------------------------------------------------
  // EVENT LOOP  
  //-----------------------------------------------------
  Bool_t eventAdded;
  if(flowControl!=fParallelCat && flowControl!=fParallelCalc) {

    ev->whichHelicityMC = whichHelicityMC;

    printf("begin loop through %lld events...\n",ev->ENT);
    for(int i=0; i<ev->ENT; i++) {

      ev->GetEvent(i);

      if(ev->Valid()) {

        // fill asymmetry plots; Asymmetry::AddEvent() checks the bin,
        // and fills plots if it's the correct bin; if not, AddEvent() does nothing
        // and returns false
        for(Int_t bn : binVec) {

          A = asymMap.at(bn);

          // set kinematic vars
          A->Mh = ev->Mh;
          A->x = ev->x;
          A->z = ev->Zpair;
          A->PhiH = ev->PhiH;
          A->PhiR = ev->PhiR;
          A->PhPerp = ev->PhPerp;
          A->Ph = ev->Ph;
          A->Q2 = ev->Q2;
          A->theta = ev->theta;

          // set spin state
          A->spinn = ev->SpinState();

          // set kinematic factors
          A->kfA = ev->GetKinematicFactor('A');
          A->kfC = ev->GetKinematicFactor('C');
          A->kfW = ev->GetKinematicFactor('W');


          // set any test modulation variables
          A->PhiTest = UNDEF;
          /*
          if(whichModulation == Asymmetry::modTest) {
            A->PhiTest = ev->GetDihadronObj()->GetSingleHadronPhiH(qA);
            A->z = ev->Z[qA];
          };
          */
          

          // add this event to Asymmetry data structures (event won't be added the bin
          // is not correct)
          eventAdded = A->AddEvent(); 

          //if(eventAdded && A->debug) ev->PrintEvent();
        };

      };
    };
  };
  // end event loop -------------------------------------------


  // concatenate spinroot files data, or just read in the concatenated file
  TSystemDirectory * sysDir;
  TSystemFile * sysFile;
  TList * sysFileList;
  TFile * appFile;
  TString appFileName;
  if(flowControl==fParallelCat) {
    sysDir = new TSystemDirectory("spinroot","spinroot");
    sysFileList = sysDir->GetListOfFiles();
    TIter nxt(sysFileList);
    while(( sysFile = (TSystemFile*) nxt() )) {
      appFileName = "spinroot/" + TString(sysFile->GetName());
      if(!sysFile->IsDirectory() && appFileName.EndsWith(".root")) {
        Tools::PrintSeparator(40,".");
        printf("concatenating data from %s\n",appFileName.Data());
        appFile = new TFile(appFileName,"READ");
        for(Int_t bn : binVec) {
          A = asymMap.at(bn);
          A->AppendData(appFile);
        };
        appFile->Close("R");
        printf("done reading %s\n",appFileName.Data());
      };
    };
    Tools::PrintSeparator(40,".");
    printf("spinroot files concatenated\n\n");
    spinrootFile->cd();
  }
  else if(flowControl==fParallelCalc) {
    appFile = new TFile("spinrootCat.root","READ");
    for(Int_t bn : binVec) {
      A = asymMap.at(bn);
      A->AppendData(appFile);
    };
  };



  // if we are filling a single spinroot file in a parallelized analysis,
  // write it out and exit
  // - this doesn't work, unfortunately. The Asymmetry objects get written, and have a
  //   data size that indicates histograms/data structures are being written, however
  //   upon trying to access anything in the TFile, it immediately seg faults
  // - alternative implementation below writes out only pertinent members of Asymmetry
  //   instead of the full Asymmetry object
  /*
  TString Aname;
  if(flowControl==fParallelFill) {
    spinrootFile->cd();
    BS->Write("BS");
    for(Int_t bn : binVec) {
      A = asymMap.at(bn);
      Aname = "A" + A->binN;
      printf("write %s\n",Aname.Data());
      A->PrintSettings();
      //A->rfData->Write(Aname);
      A->Write(Aname);
    };
    spinrootFile->Close();
    return 1;
  };
  */
  if(flowControl==fParallelFill || flowControl==fParallelCat) {
    spinrootFile->cd();
    BS->Write("BS");
    for(Int_t bn : binVec) {
      A = asymMap.at(bn);
      A->StreamData(spinrootFile);
    };
    spinrootFile->Close();
    return 1;
  };
    


  // compute asymmetries
  Float_t asymValue,asymError;
  Float_t RFasymValue[Asymmetry::nAmp];
  Float_t RFasymError[Asymmetry::nAmp];
  Float_t meanKF;
  Float_t kinValue,kinError;
  Float_t chisq,ndf;
  printf("--- calculate asymmetries\n");
  for(Int_t bn : binVec) {
    A = asymMap.at(bn);
    A->CalculateAsymmetries();

    if( ( A->asym2d==false && A->fitFunc!=NULL) || 
        ( A->asym2d==true && A->fitFunc2!=NULL) ) {

      // asymmetry value and statistical uncertainty
      if(!(A->asym2d)) {
        asymValue = A->fitFunc->GetParameter(1);
        asymError = A->fitFunc->GetParError(1);
      } else {
        asymValue = A->fitFunc2->GetParameter(0);
        asymError = A->fitFunc2->GetParError(0);
      };
      for(int aa=0; aa<N_AMP; aa++) {
        RFasymValue[aa] = A->rfA[aa]->getVal();
        RFasymError[aa] = A->rfA[aa]->getError();
      };


      // divide out mean kinematic factor
      /*
      meanKF = A->kfDist->GetMean();
      asymValue /= meanKF;
      for(int aa=0; aa<N_AMP; aa++) RFasymValue[aa] /= meanKF;
      */



      // IV value and uncertainty
      switch(dimensions) {
        case 1:
          kinValue = A->ivDist1->GetMean();
          kinError = A->ivDist1->GetRMS();
          break;
        case 2:
          kinValue = A->ivDist2->GetMean(1);
          kinError = A->ivDist2->GetRMS(1);
          break;
        case 3:
          kinValue = A->ivDist3->GetMean(1);
          kinError = A->ivDist3->GetRMS(1);
          break;
      };
      kinError = 0; // OVERRIDE

      // chi2 and ndf
      if(!(A->asym2d)) {
        chisq = A->fitFunc->GetChisquare();
        ndf = A->fitFunc->GetNDF();
      } else {
        chisq = A->fitFunc2->GetChisquare();
        ndf = A->fitFunc2->GetNDF();
      };
      chisqDist->Fill(chisq);

      // set points
      kindepGr = kindepMap.at(bn);
      kindepGr->SetPoint(A->B[0],kinValue,asymValue);
      kindepGr->SetPointError(A->B[0],kinError,asymError);

      for(int aa=0; aa<N_AMP; aa++) {
        RFkindepGr[aa] = RFkindepMap[aa].at(bn);
        RFkindepGr[aa]->SetPoint(A->B[0],kinValue,RFasymValue[aa]);
        RFkindepGr[aa]->SetPointError(A->B[0],kinError,RFasymError[aa]);
      };

      chindfGr = chindfMap.at(bn);
      chindfGr->SetPoint(A->B[0],kinValue,chisq/ndf);

      rellumGr = rellumMap.at(bn);
      rellumGr->SetPoint(A->B[0],kinValue,A->rellum);
      rellumGr->SetPointError(A->B[0],0,A->rellumErr);

    };

  };



  // -- instantiate canvases
  TString canvNameSuffix = "Canv_" + modN;
  TString canvName;
  for(int d=0; d<dimensions; d++) {
    if(d==1) canvNameSuffix += "_bins_" + BS->IVname[ivVar[d]];
    else canvNameSuffix += "_" + BS->IVname[ivVar[d]];
  };
  Int_t canvX,canvY,divX,divY;
  Int_t canvModX,canvModY,divModX,divModY;
  Int_t canvSize = 800;
  switch(dimensions) {
    case 1:
      canvX=canvSize; canvY=canvSize;
      divX=1; divY=1;
      canvModX=NB[0]*canvSize; canvModY=canvSize;
      divModX=NB[0]; divModY=1;
      break;
    case 2:
      canvX=NB[1]*canvSize; canvY=canvSize;
      divX=NB[1]; divY=1;
      canvModX=NB[1]*canvSize; canvModY=NB[0]*canvSize;
      divModX=NB[1]; divModY=NB[0];
      break;
    case 3:
      canvX=NB[2]*canvSize; canvY=NB[1]*canvSize;
      divX=NB[2]; divY=NB[1];
      canvModX=canvSize; canvModY=canvSize; // (not used) 
      divModX=1; divModY=1; // (not used) 
      break;
  };

  canvName = "kindep" + canvNameSuffix;
  TCanvas * kindepCanv = new TCanvas(canvName,canvName,canvX,canvY);
  kindepCanv->Divide(divX,divY);

  TCanvas * RFkindepCanv[Asymmetry::nAmp];
  for(int aa=0; aa<N_AMP; aa++) {
    canvName = "RF_A"+TString::Itoa(aa,10)+"_kindep" + canvNameSuffix;
    RFkindepCanv[aa] = new TCanvas(canvName,canvName,canvX,canvY);
    RFkindepCanv[aa]->Divide(divX,divY);
  };

  canvName = "chindf" + canvNameSuffix;
  TCanvas * chindfCanv = new TCanvas(canvName,canvName,canvX,canvY);
  chindfCanv->Divide(divX,divY);

  canvName = "rellum" + canvNameSuffix;
  TCanvas * rellumCanv = new TCanvas(canvName,canvName,canvX,canvY);
  rellumCanv->Divide(divX,divY);

  canvName = "asymMod" + canvNameSuffix;
  TCanvas * asymModCanv = new TCanvas(canvName,canvName,canvModX,canvModY); 
  asymModCanv->Divide(divModX,divModY);

  canvName = "asymModHist2" + canvNameSuffix;
  TCanvas * asymModHist2Canv = new TCanvas(canvName,canvName,canvModX,canvModY); 
  asymModHist2Canv->Divide(divModX,divModY);

  canvName = "modDist" + canvNameSuffix;
  TCanvas * modDistCanv = new TCanvas(canvName,canvName,canvModX,canvModY); 
  modDistCanv->Divide(divModX,divModY);




  // -- add objects to canvases and graphs to multiGr
  Int_t pad;
  TGraphErrors * RFkindepGrClone[Asymmetry::nAmp];
  if(dimensions==1) {
    binNum = HashBinNum(0);

    kindepGr = kindepMap.at(binNum);
    kindepCanv->cd();
    DrawKinDepGraph(kindepGr,BS,ivVar[0]);

    multiGr = multiMap.at(binNum);
    multiGr->Add(kindepGr); // include linear fit in multiGr
    for(int aa=0; aa<N_AMP; aa++) {
      RFkindepGr[aa] = RFkindepMap[aa].at(binNum);
      RFkindepCanv[aa]->cd();
      DrawKinDepGraph(RFkindepGr[aa],BS,ivVar[0]);
      RFkindepGrClone[aa] = ShiftGraph(RFkindepGr[aa],aa+1);
      multiGr->Add(RFkindepGrClone[aa]);
    };

    chindfGr = chindfMap.at(binNum);
    chindfCanv->cd();
    DrawSimpleGraph(chindfGr,BS,ivVar[0]);

    rellumGr = rellumMap.at(binNum);
    rellumCanv->cd();
    DrawSimpleGraph(rellumGr,BS,ivVar[0],false);

    for(int b0=0; b0<NB[0]; b0++) {
      binNum = HashBinNum(b0);
      A = asymMap.at(binNum);
      asymModCanv->cd(b0+1);
      if(!(A->asym2d)) DrawAsymGr(A->asymGr);
      else DrawAsymGr2(A->asymGr2);
      modDistCanv->cd(b0+1);
      if(!(A->asym2d)) A->modDist->Draw();
      else A->modDist2->Draw("colz");
      if(A->asym2d) {
        asymModHist2Canv->cd(b0+1);
        A->asymGr2hist->Draw("colz");
      };
    };
  }
  else if(dimensions==2) {
    for(int b1=0; b1<NB[1]; b1++) {
      pad = b1+1;
      binNum = HashBinNum(0,b1);

      kindepGr = kindepMap.at(binNum);
      kindepCanv->cd(pad);
      DrawKinDepGraph(kindepGr,BS,ivVar[0]);

      multiGr = multiMap.at(binNum);
      multiGr->Add(kindepGr); // include linear fit in multiGr

      for(int aa=0; aa<N_AMP; aa++) {
        RFkindepGr[aa] = RFkindepMap[aa].at(binNum);
        RFkindepCanv[aa]->cd(pad);
        DrawKinDepGraph(RFkindepGr[aa],BS,ivVar[0]);
        RFkindepGrClone[aa] = ShiftGraph(RFkindepGr[aa],aa+1);
        multiGr->Add(RFkindepGrClone[aa]);
      };

      chindfGr = chindfMap.at(binNum);
      chindfCanv->cd(pad);
      DrawSimpleGraph(chindfGr,BS,ivVar[0]);

      rellumGr = rellumMap.at(binNum);
      rellumCanv->cd(pad);
      DrawSimpleGraph(rellumGr,BS,ivVar[0],false);

      for(int b0=0; b0<NB[0]; b0++) {
        binNum = HashBinNum(b0,b1);
        A = asymMap.at(binNum);
        asymModCanv->cd(b0*NB[1]+b1+1);
        if(!(A->asym2d)) DrawAsymGr(A->asymGr);
        else DrawAsymGr2(A->asymGr2);
        modDistCanv->cd(b0*NB[1]+b1+1);
        if(!(A->asym2d)) A->modDist->Draw();
        else A->modDist2->Draw();
        if(A->asym2d) {
          asymModHist2Canv->cd(b0*NB[1]+b1+1);
          A->asymGr2hist->Draw("colz");
        };
      };
    };
  }
  else if(dimensions==3) {
    for(int b1=0; b1<NB[1]; b1++) {
      for(int b2=0; b2<NB[2]; b2++) {
        pad = b1*NB[2]+b2+1;
        binNum = HashBinNum(0,b1,b2);

        kindepGr = kindepMap.at(binNum);
        kindepCanv->cd(pad);
        DrawKinDepGraph(kindepGr,BS,ivVar[0]);

        multiGr = multiMap.at(binNum);
        multiGr->Add(kindepGr); // include linear fit in multiGr

        for(int aa=0; aa<N_AMP; aa++) {
          RFkindepGr[aa] = RFkindepMap[aa].at(binNum);
          RFkindepCanv[aa]->cd(pad);
          DrawKinDepGraph(RFkindepGr[aa],BS,ivVar[0]);
          RFkindepGrClone[aa] = ShiftGraph(RFkindepGr[aa],aa+1);
          multiGr->Add(RFkindepGrClone[aa]);
        };

        chindfGr = chindfMap.at(binNum);
        chindfCanv->cd(pad);
        DrawSimpleGraph(chindfGr,BS,ivVar[0]);

        rellumGr = rellumMap.at(binNum);
        rellumCanv->cd(pad);
        DrawSimpleGraph(rellumGr,BS,ivVar[0],false);
      };
    };
  };


  // sum distributions (for showing bin boundaries)
  TH1D * ivFullDist1;
  TH2D * ivFullDist2;
  TH3D * ivFullDist3;
  TH1D * modFullDist;
  TH2D * modFullDist2; // for 2d modulation
  TH2D * IVvsModFullDist;

  if(dimensions==1) {
    for(int b0=0; b0<NB[0]; b0++) {
      binNum = HashBinNum(b0);
      A = asymMap.at(binNum);
      if(b0==0) {
        ivFullDist1 = (TH1D*)(A->ivDist1)->Clone();
        SetCloneName(ivFullDist1);
        if(!(A->asym2d)) {
          IVvsModFullDist = (TH2D*)(A->IVvsModDist)->Clone();
          SetCloneName(IVvsModFullDist);
          modFullDist = (TH1D*)(A->modDist)->Clone();
          SetCloneName(modFullDist);
        } else {
          modFullDist2 = (TH2D*)(A->modDist2)->Clone();
          SetCloneName(modFullDist2);
        };
      } else {
        ivFullDist1->Add(A->ivDist1);
        if(!(A->asym2d)) {
          IVvsModFullDist->Add(A->IVvsModDist);
          modFullDist->Add(A->modDist);
        } else {
          modFullDist2->Add(A->modDist2);
        };
      };
    };
  }
  else if(dimensions==2) {
    for(int b1=0; b1<NB[1]; b1++) {
      for(int b0=0; b0<NB[0]; b0++) {
        binNum = HashBinNum(b0,b1);
        A = asymMap.at(binNum);
        if(b0==0 && b1==0) {
          ivFullDist2 = (TH2D*)(A->ivDist2)->Clone();
          SetCloneName(ivFullDist2);
          if(!(A->asym2d)) {
            modFullDist = (TH1D*)(A->modDist)->Clone();
            SetCloneName(modFullDist);
          } else {
            modFullDist2 = (TH2D*)(A->modDist2)->Clone();
            SetCloneName(modFullDist2);
          };
        } else {
          ivFullDist2->Add(A->ivDist2);
          if(!(A->asym2d)) modFullDist->Add(A->modDist);
          else modFullDist2->Add(A->modDist2);
        };
      };
    };
  }
  else if(dimensions==3) {
    for(int b2=0; b2<NB[2]; b2++) {
      for(int b1=0; b1<NB[1]; b1++) {
        for(int b0=0; b0<NB[0]; b0++) {
          binNum = HashBinNum(b0,b1,b2);
          A = asymMap.at(binNum);
          if(b0==0 && b1==0 && b2==0) {
            ivFullDist3 = (TH3D*)(A->ivDist3)->Clone();
            SetCloneName(ivFullDist3);
            if(!(A->asym2d)) {
              modFullDist = (TH1D*)(A->modDist)->Clone();
              SetCloneName(modFullDist);
            } else {
              modFullDist2 = (TH2D*)(A->modDist2)->Clone();
              SetCloneName(modFullDist2);
            };
          } else {
            ivFullDist3->Add(A->ivDist3);
            if(!(A->asym2d)) modFullDist->Add(A->modDist);
            else modFullDist2->Add(A->modDist2);
          };
        };
      };
    };
  };



  // write output to TFile
  if(flowControl==fSerial || 
     flowControl==fSerialRenamed || 
     flowControl==fParallelCalc
  ) {

    resultFile->cd();
    // -- Asymmetry objects
    printf("--- write Asymmetry objects\n");
    for(Int_t bn : binVec) {
      A = asymMap.at(bn);
      A->StreamData(resultFile);
    };

    // -- "full" distributions
    if(dimensions==1) {
      ivFullDist1->Write();
      if(!(A->asym2d)) IVvsModFullDist->Write();
    } else if(dimensions==2) {
      ivFullDist2->Write();
    } else if(dimensions==3) {
      ivFullDist3->Write();
    };
    if(!(A->asym2d)) modFullDist->Write();
    else modFullDist2->Write();


    // -- asymmetries and kindep graphs
    printf("--- write kinematic-dependent asymmetries\n");
    for(Int_t bn : binVec) {
      A = asymMap.at(bn);
      A->PrintSettings();

      // first write out the asymmetry vs. modulation graphs
      if(!(A->asym2d)) A->asymGr->Write();
      else A->asymGr2->Write();

      // then write out the kindepGr *after* writing out all the
      // relevant asymmetry vs. modulation graphs
      if(A->B[0] + 1 == NB[0]) {
        binNum = HashBinNum(A->B[0], A->B[1], A->B[2]);
        kindepGr = kindepMap.at(binNum);
        kindepGr->Write();
        for(int aa=0; aa<N_AMP; aa++) {
          RFkindepGr[aa] = RFkindepMap[aa].at(binNum);
          RFkindepGr[aa]->Write();
        };
        multiGr = multiMap.at(binNum);
        multiGr->Write();
      };
    };

    kindepCanv->Write();
    for(int aa=0; aa<N_AMP; aa++) RFkindepCanv[aa]->Write();
    chindfCanv->Write();
    chisqDist->Write();
    rellumCanv->Write();
    if(dimensions==1 || dimensions==2) {
      asymModCanv->Write();
      if(A->asym2d) asymModHist2Canv->Write();
      modDistCanv->Write();
    };


    // -- RooFit results
    TCanvas * rfCanv[Asymmetry::nAmp];
    TString rfCanvName[Asymmetry::nAmp];

    for(Int_t bn : binVec) {
      A = asymMap.at(bn);

      for(int aa=0; aa<N_AMP; aa++) {
        rfCanvName[aa] = "RF_A" + TString::Itoa(aa,10) + "_NLL_" + A->binN;
        rfCanv[aa] = new TCanvas(rfCanvName[aa],rfCanvName[aa],800,800);
        A->rfNLLplot[aa]->Draw();
        rfCanv[aa]->Write();
      };

      Tools::PrintTitleBox("roofit function");
      A->PrintSettings();
      printf("\n");
      for(int ss=0; ss<nSpin; ss++) {
        printf("%s: %s\n",SpinTitle(ss).Data(),A->rfPdfFormu[ss].Data());
      };
      printf("\n");

      printf("fit parameter results:\n");
      for(int aa=0; aa<N_AMP; aa++) {
        printf(" >> A%d = %.3f +/- %.3f\n",
            aa, A->rfA[aa]->getVal(), A->rfA[aa]->getError() );
      };
      for(int dd=0; dd<N_D; dd++) {
        printf(" >> D%d = %.3f +/- %.3f\n",
            dd, A->rfD[dd]->getVal(), A->rfD[dd]->getError() );
      };
      /*
         printf(" >> Y+ = %.3f +/- %.3f\n",
         A->rfYield[0]->getVal(), A->rfYield[0]->getError() );
         printf(" >> Y- = %.3f +/- %.3f\n",
         A->rfYield[1]->getVal(), A->rfYield[1]->getError() );
         */

      printf("\n");

    };


    // DEPRECATED
    // print modDist boundaries (used for determining modDist boundaries
    // which depend on kinematics, for g1perp modulation PhPerp/Mh-scaling test)
    /*
       Float_t mbound;
       if(dimensions==1 && whichModulation==Asymmetry::scaleSinPhiHR) {
       printf("MBOUND\n");
       printf("if(v_==v%s) {\n",(BS->IVname[ivVar[0]]).Data());
       for(int b=0; b<NB[0]; b++) {
       binNum = HashBinNum(b);
       A = asymMap.at(binNum);
       mbound = TMath::Max(
       TMath::Abs(Tools::GetFirstFilledX(A->modDist)),
       TMath::Abs(Tools::GetLastFilledX(A->modDist))
       );
       mbound += 0.1;
       printf("  if(b_==%d) return %f;\n",b,mbound);
       };
       printf("};\n");
       };
       */



    // print images
    TString pngName;
    if(flowControl==fSerialRenamed) {
      pngName = Form("spinout/%s.%s.png",kindepCanv->GetName(),dihName.Data());
      kindepCanv->Print(pngName,"png");
      for(int aa=0; aa<N_AMP; aa++) {
        pngName = Form("spinout/%s.%s.png",RFkindepCanv[aa]->GetName(),dihName.Data());
        RFkindepCanv[aa]->Print(pngName,"png");
      };
      pngName = Form("spinout/%s.%s.png",chindfCanv->GetName(),dihName.Data());
      chindfCanv->Print(pngName,"png");
      pngName = Form("spinout/%s.%s.png",rellumCanv->GetName(),dihName.Data());
      rellumCanv->Print(pngName,"png");
      if(dimensions==1 || dimensions==2) {
        pngName = Form("spinout/%s.%s.png",asymModCanv->GetName(),dihName.Data());
        asymModCanv->Print(pngName,"png"); 
        pngName = Form("spinout/%s.%s.png",modDistCanv->GetName(),dihName.Data());
        modDistCanv->Print(pngName,"png"); 
      };
    };

  };



  if(flowControl==fSerial || 
     flowControl==fSerialRenamed || 
     flowControl==fParallelCalc
  ) {
    resultFile->Close();
  };

  printf("--- end %s\n",argv[0]);

  return 0;
};


Int_t HashBinNum(Int_t bin0, Int_t bin1, Int_t bin2) {
  Int_t retval = bin0;
  if(bin1>=0) retval += 10 * bin1;
  if(bin2>=0) retval += 100 * bin2;
  return retval;
};


void DrawKinDepGraph(TGraphErrors * g_, Binning * B_, Int_t v_) {

  g_->Draw("APE"); // draw once, so we can then format it

  g_->SetLineWidth(2);

  g_->SetMarkerStyle(kFullCircle);
  g_->SetMarkerColor(kBlack);
  g_->SetMarkerSize(1.3);

  // set vertical axis range (it is overridden if the plot's vertical range
  // is larger than the desired range)
  Float_t yMin = -0.12;
  Float_t yMax = 0.12;
  if(g_->GetYaxis()->GetXmin() < yMin) yMin = g_->GetYaxis()->GetXmin();
  if(g_->GetYaxis()->GetXmax() > yMax) yMax = g_->GetYaxis()->GetXmax();
  g_->GetYaxis()->SetRangeUser(yMin,yMax);

  // set horizontal range
  //g_->GetXaxis()->SetLimits(B_->minIV[v_],B_->maxIV[v_]);

  g_->Draw("APE"); // draw again to apply the formatting


  // zero line
  Float_t drawMin = g_->GetXaxis()->GetXmin();
  Float_t drawMax = g_->GetXaxis()->GetXmax();
  TLine * zeroLine = new TLine(drawMin,0,drawMax,0);
  zeroLine->SetLineColor(kBlack);
  zeroLine->SetLineWidth(1.5);
  zeroLine->SetLineStyle(kDashed);
  zeroLine->Draw();
};


void DrawSimpleGraph(TGraphErrors * g_, Binning * B_, Int_t v_, Bool_t setRange) {

  g_->Draw("AP"); // draw once, so we can then format it

  //g_->SetLineWidth(2);

  g_->SetMarkerStyle(kFullCircle);
  g_->SetMarkerColor(kBlack);
  g_->SetMarkerSize(1.3);

  if(setRange) {
    // set vertical axis range (it is overridden if the plot's vertical range
    // is larger than the desired range)
    Float_t yMin = 0;
    Float_t yMax = 2;
    if(g_->GetYaxis()->GetXmin() < yMin) yMin = g_->GetYaxis()->GetXmin() - 0.2;
    if(g_->GetYaxis()->GetXmax() > yMax) yMax = g_->GetYaxis()->GetXmax() + 0.2;
    g_->GetYaxis()->SetRangeUser(yMin,yMax);

    // set horizontal range
    //g_->GetXaxis()->SetLimits(B_->minIV[v_],B_->maxIV[v_]);
  };

  g_->Draw("AP"); // draw again to apply the formatting


  // unity line
  Float_t drawMin = g_->GetXaxis()->GetXmin();
  Float_t drawMax = g_->GetXaxis()->GetXmax();
  TLine * unityLine = new TLine(drawMin,1,drawMax,1);
  unityLine->SetLineColor(kBlack);
  unityLine->SetLineWidth(1.5);
  unityLine->SetLineStyle(kDashed);
  unityLine->Draw();
};


void DrawAsymGr(TGraphErrors * g_) {

  TString titleTmp = g_->GetTitle();
  g_->SetTitle(TString(dihTitle+" "+titleTmp));

  g_->Draw("APE"); // draw once, so we can then format it

  g_->SetLineColor(kBlack);
  g_->SetLineWidth(2);

  g_->SetMarkerStyle(kFullCircle);
  g_->SetMarkerColor(kBlack);
  g_->SetMarkerSize(1.3);

  // set vertical axis range (it is overridden if the plot's vertical range
  // is larger than the desired range)
  Float_t yMin = -0.2;
  Float_t yMax = 0.2;
  if(g_->GetYaxis()->GetXmin() < yMin) yMin = g_->GetYaxis()->GetXmin() - 0.05;
  if(g_->GetYaxis()->GetXmax() > yMax) yMax = g_->GetYaxis()->GetXmax() + 0.05;
  g_->GetYaxis()->SetRangeUser(yMin,yMax);

  g_->Draw("APE"); // draw again to apply the formatting

};


void DrawAsymGr2(TGraph2DErrors * g_) {

  TString titleTmp = g_->GetTitle();
  g_->SetTitle(TString(dihTitle+" "+titleTmp));

  g_->Draw("ERR P"); // draw once, so we can then format it

  g_->SetLineColor(kBlack);
  g_->SetLineWidth(2);

  g_->SetMarkerStyle(kFullCircle);
  g_->SetMarkerColor(kBlack);
  g_->SetMarkerSize(1.3);

  // set vertical axis range (it is overridden if the plot's vertical range
  // is larger than the desired range)
  /*
     Float_t yMin = -0.2;
     Float_t yMax = 0.2;
     if(g_->GetYaxis()->GetXmin() < yMin) yMin = g_->GetYaxis()->GetXmin() - 0.05;
     if(g_->GetYaxis()->GetXmax() > yMax) yMax = g_->GetYaxis()->GetXmax() + 0.05;
     g_->GetYaxis()->SetRangeUser(yMin,yMax);
     */

  g_->Draw("ERR P"); // draw again to apply the formatting

};


void SetCloneName(TH1 * clone_) {
  TString cloneName,cloneTitle;
  cloneName = clone_->GetName();
  cloneName.ReplaceAll("Dist","FullDist");
  cloneName.ReplaceAll("0","");
  cloneTitle = clone_->GetTitle();
  cloneTitle(TRegexp("::.*$")) = "";
  clone_->SetName(cloneName);
  clone_->SetTitle(cloneTitle);
};


// shift a graph's points to the right slightly (a clone of the graph, with shifted
// points, is returned)
TGraphErrors * ShiftGraph(TGraphErrors * gr, Int_t nShift) {
  //TGraphErrors * retGr = (TGrapErrors*) gr->Clone();
  TGraphErrors * retGr = new TGraphErrors();
  Double_t * grX = gr->GetX();
  Double_t * grY = gr->GetY();
  Double_t * grEX = gr->GetEX();
  Double_t * grEY = gr->GetEY();
  for(int nn=0; nn<gr->GetN(); nn++) {
    //retGr->SetPoint(nn, grX[nn]+nShift*0.01, grY[nn]); // shifting enabled
    retGr->SetPoint(nn, grX[nn], grY[nn]); // shifting disabled
    retGr->SetPointError(nn, grEX[nn], grEY[nn]);
  };

  switch(nShift) {
    case 1:
      retGr->SetLineColor(N_AMP==1?kGray+3:kGreen+1); 
      retGr->SetLineStyle(2);
      retGr->SetMarkerStyle(22);
      break;
    case 2:
      retGr->SetLineColor(kRed); 
      retGr->SetLineStyle(1);
      retGr->SetMarkerStyle(kFullCircle);
      break;
    case 3:
      retGr->SetLineColor(kBlue);
      retGr->SetLineStyle(3);
      retGr->SetMarkerStyle(23);
      break;
    case 4:
      retGr->SetLineColor(kMagenta);
      retGr->SetMarkerStyle(kFullCircle);
      break;
    default: retGr->SetLineColor(kGray);
  };
  

  retGr->SetMarkerColor(kBlack);
  retGr->SetLineWidth(2);
  retGr->SetMarkerSize(1.3);

  return retGr;
};


// set default arguments
void SetDefaultArgs() {
  inputData = "";
  pairType = EncodePairType(kPip,kPim);
  whichModulation = Asymmetry::weightSinPhiHR;
  ivType = Binning::vM + 1;
  flowControl = fSerial;
  whichHelicityMC = 0;

  DecodePairType(pairType,whichHad[qA],whichHad[qB]);
  dihTitle = PairTitle(whichHad[qA],whichHad[qB]);
};


// help printout
int PrintUsage() {

  SetDefaultArgs();
  BS = new Binning(EncodePairType(kPip,kPim));
  fprintf(stderr,"\nUSAGE: asym.exe [-f or -d input_data ] [options...]\n\n");

  printf("INPUT DATA:\n");
  printf(" -f\tsingle ROOT file\n");
  printf(" -d\tdirectory of ROOT files\n");
  printf(" NOTE: specify input with either -f or -d, but not both\n");
  printf("\n");

  printf("OPTIONS:\n");

  printf(" -p\tpair type, specified as a hexadecimal number\n");
  printf("   \trun PrintEnumerators.C for notation\n");
  printf("   \tdefault = 0x%x (%s)\n\n",pairType,dihTitle.Data());

  printf(" -m\tazimuthal modulation for asymmetry\n");
  for(int m=0; m<Asymmetry::nMod; m++) {
    A = new Asymmetry(BS,m,UNDEF);
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

  printf(" -c\tflow control (typically used by wrapper scripts)\n");
  printf("   \t 0 = default behavior\n");
  printf("   \t 1 = rename output spin.root rootfile and prints pngs\n");
  printf("   \t 2,3 = for parallel processing\n");
  printf("   \tdefault = %d\n\n",flowControl);

  printf(" -g\t (for MC) - select which helicityMC to use\n\n");

  return 0;
};
