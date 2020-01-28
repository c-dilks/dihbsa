// asymmetry fitter
// - fits spinroot cat file (from catSpinroot.cpp) for extracting asymmetry amplitudes
// - produces the file spinroot/asym.root
// - you can specify a specific `spinroot` directory, if you want

#include <cstdlib>
#include <iostream>
#include <vector>
#include <map>

// ROOT
#include "TFile.h"
#include "TString.h"
#include "TMath.h"
#include "TRegexp.h"
#include "TH1.h"
#include "TH2.h"
#include "TF1.h"
#include "TCanvas.h"
#include "TLine.h"
#include "TStyle.h"
#include "TMultiGraph.h"
#include "TSystem.h"
#include "TObjArray.h"
#include "TLegend.h"

// DihBsa
#include "Constants.h"
#include "Binning.h"
#include "Asymmetry.h"
#include "Modulation.h"


// subroutines
void DrawKinDepGraph(TGraphErrors * g_, Binning * B_, Int_t d_);
void DrawSimpleGraph(TGraphErrors * g_, Binning * B_, Int_t d_, Bool_t setRange_=true);
void DrawAsymGr(TGraphErrors * g_);
void DrawAsymGr2(TGraph2DErrors * g_);
TGraphErrors * ShiftGraph(TGraphErrors * gr, Int_t nShift);
void SetCloneName(TH1 * clone_);

// global variables
Int_t N_AMP,N_D;
TString dihTitle,dihName;
Binning * BS;
Asymmetry * A;
Modulation * modu;



int main(int argc, char** argv) {

  //////////////////////////////////////////////
  // ARGUMENTS
  TString spinrootDir = "spinroot";
  Int_t fitMode = 0;
  Float_t DparamVal = 0; // (for systematic uncertainty from D_1 pp-wave)
  if(argc>1) fitMode = (Int_t)strtof(argv[1],NULL);
  if(argc>2) spinrootDir = TString(argv[2]);
  if(argc>3) DparamVal = (Float_t)strtof(argv[3],NULL);
  //////////////////////////////////////////////


  //////////////////////////////////////////////
  // OPTIONS
  Bool_t includeOAonMultiGr = false;
  gStyle->SetOptFit(1);
  //////////////////////////////////////////////


  // open spinroot cat file and result file
  TString asymFileN = Form("%s/asym_%d.root",spinrootDir.Data(),fitMode);
  TFile * asymFile = new TFile(asymFileN,"RECREATE");
  TFile * catFile = new TFile(TString(spinrootDir+"/cat.root"),"READ");


  // instantiate Binning and Asymmetry, and extract them from catFile
  std::map<Int_t, Asymmetry*> asymMap;
  catFile->GetObject("BS",BS);
  for(Int_t bn : BS->binVec) {
    A = new Asymmetry(BS,bn);
    if(A->success) {
      A->AppendData(catFile);
      asymMap.insert(std::pair<Int_t,Asymmetry*>(bn,A));
    }
    else return 0;
  };
  dihTitle = PairTitle(BS->whichHad[qA],BS->whichHad[qB]);
  dihName = PairName(BS->whichHad[qA],BS->whichHad[qB]);
  modu = new Modulation();


  // print which IV will be analyzed
  TString modN = A->oaModulationName;
  printf("--------> Analysing %s asymmetries vs. %s ",
    modN.Data(),(BS->GetIVname(0)).Data());
  if(BS->dimensions>=2) printf("in bins of %s ",(BS->GetIVname(1)).Data());
  if(BS->dimensions>=3) printf("and %s ",(BS->GetIVname(2)).Data());
  printf("\n\n");


  // perform the fits
  printf("--- calculate asymmetries\n");
  for(Int_t bn : BS->binVec) {
    A = asymMap.at(bn);
    A->FitOneAmp();
    A->FitMultiAmp(fitMode,DparamVal);
  };


  // build maps from Binning::binVec number to plots etc.
  TGraphErrors * kindepGr;
  TGraphErrors * RFkindepGr[Asymmetry::nAmp];
  TGraphErrors * chindfGr;
  TGraphErrors * rellumGr;
  TMultiGraph * multiGr;

  std::map<Int_t, TGraphErrors*> kindepMap;
  std::map<Int_t, TGraphErrors*> RFkindepMap[Asymmetry::nAmp];
  std::map<Int_t, TMultiGraph*> multiMap;
  std::map<Int_t, TGraphErrors*> chindfMap;
  std::map<Int_t, TGraphErrors*> rellumMap;

  TString grTitle,grTitleSuffix,grName,grNameSuffix;
  Bool_t first = true;
  for(Int_t bn : BS->binVec) {
    A = asymMap.at(bn);

    // get number of RooFit params
    if(first) { 
      N_AMP = A->nAmpUsed;
      N_D = A->nDparamUsed;
      first = false;
    };

    // set graph title and name suffixes
    switch(BS->dimensions) {
      case 1:
        grTitleSuffix = BS->GetIVtitle(0) + ";" + BS->GetIVtitle(0);
        grNameSuffix = BS->GetIVname(0);
        break;
      case 2:
        grTitleSuffix = BS->GetIVtitle(0) + " :: " +
          BS->GetBoundStr(bn,1) + ";" +
          BS->GetIVtitle(0);
        grNameSuffix = Form("%s_bin_%s%d",
            (BS->GetIVname(0)).Data(),
            (BS->GetIVname(1)).Data(), BS->UnhashBinNum(bn,1)
            );
        break;
      case 3:
        grTitleSuffix = BS->GetIVtitle(0) + " :: " +
          BS->GetBoundStr(bn,1) + ", " +
          BS->GetBoundStr(bn,2) + ";" +
          BS->GetIVtitle(0);
        grNameSuffix = Form("%s_bin_%s%d_%s%d",
            (BS->GetIVname(0)).Data(),
            (BS->GetIVname(1)).Data(), BS->UnhashBinNum(bn,1),
            (BS->GetIVname(2)).Data(), BS->UnhashBinNum(bn,2)
            );
        break;
    };


    // instantiate graphs; only needs to be done for each IV1 and IV2 bin, since the
    // horizontal axis of these graphs are all IV0 (hence the if statement here)
    if(BS->UnhashBinNum(bn,0)==0) {

      // instantiate kindep graph, for linear fit ("l.f.") result
      grTitle = dihTitle + " A_{LU}[" + A->oaModulationTitle + "] " + 
        " vs. " + grTitleSuffix;
      grName = "kindep_" + grNameSuffix;
      kindepGr = new TGraphErrors();
      kindepGr->SetName(grName);
      kindepGr->SetTitle(grTitle);

      // instantiate multiGraph, for plotting kindep graphs together
      grTitle = dihTitle + " A_{LU}[" + A->oaModulationTitle + "] " + 
        " vs. " + grTitleSuffix;
      grName = "multiGr_" + grNameSuffix;
      multiGr = new TMultiGraph();
      multiGr->SetName(grName);
      multiGr->SetTitle(grTitle);

      // instantiate kindep graphs for maximum likelihood method result
      // for each fit parameter
      for(int aa=0; aa<N_AMP; aa++) {
        grTitle = dihTitle + " " + TString(A->rfA[aa]->GetTitle()) + 
          " vs. " + grTitleSuffix;
        grName = "RF_A" + TString::Itoa(aa,10) + "_kindep_" + grNameSuffix;
        RFkindepGr[aa] = new TGraphErrors();
        RFkindepGr[aa]->SetName(grName);
        RFkindepGr[aa]->SetTitle(grTitle);
      };

      // instantiate chi2/ndf graphs
      grTitle = "#chi^{2}/NDF of " +
        dihTitle + " A_{LU}[" + A->oaModulationTitle + "]_{l.f.} " + 
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

    }


    // insert objects into maps (note: these are many-to-one maps, i.e., several
    // bin numbers will map to the same pointer)
    kindepMap.insert(std::pair<Int_t,TGraphErrors*>(bn,kindepGr));
    multiMap.insert(std::pair<Int_t,TMultiGraph*>(bn,multiGr));
    chindfMap.insert(std::pair<Int_t,TGraphErrors*>(bn,chindfGr));
    rellumMap.insert(std::pair<Int_t,TGraphErrors*>(bn,rellumGr));
    for(int aa=0; aa<N_AMP; aa++) {
      RFkindepMap[aa].insert(std::pair<Int_t,TGraphErrors*>(bn,RFkindepGr[aa]) );
    };
  };

  // overall summary plots
  TH1F * chisqDist = new TH1F("chisqDist",
      "#chi^{2} distribution (from linear fit results)",100,0,20);

  // graph asymmetry results
  Float_t asymValue,asymError;
  Float_t RFasymValue[Asymmetry::nAmp];
  Float_t RFasymError[Asymmetry::nAmp];
  Float_t meanKF;
  Float_t kinValue,kinError;
  Float_t chisq,ndf;
  for(Int_t bn : BS->binVec) {
    A = asymMap.at(bn);

    if( ( A->oa2d==false && A->fitFunc!=NULL) || 
        ( A->oa2d==true && A->fitFunc2!=NULL) ) {

      // linear fit result's asymmetry value and statistical uncertainty
      if(!(A->oa2d)) {
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
      switch(BS->dimensions) {
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
      kinError = 0; // OVERRIDE (since this should be a systematic uncertainty)

      // chi2 and ndf
      if(!(A->oa2d)) {
        chisq = A->fitFunc->GetChisquare();
        ndf = A->fitFunc->GetNDF();
      } else {
        chisq = A->fitFunc2->GetChisquare();
        ndf = A->fitFunc2->GetNDF();
      };
      chisqDist->Fill(chisq);

      // set graph points
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
  Int_t NB[3]; // # of bins for each IV
  for(int d=0; d<BS->dimensions; d++) NB[d] = BS->GetNbins(d);
  TString canvNameSuffix = "Canv_" + modN;
  TString canvName;
  for(int d=0; d<BS->dimensions; d++) {
    if(d==1) canvNameSuffix += "_bins_" + BS->GetIVname(d);
    else canvNameSuffix += "_" + BS->GetIVname(d);
  };
  Int_t canvX,canvY,divX,divY;
  Int_t canvModX,canvModY,divModX,divModY;
  Int_t canvSize = 800;
  switch(BS->dimensions) {
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

  // -- add objects to canvases and graphs to multigraphs
  Int_t pad;
  TGraphErrors * RFkindepGrClone[Asymmetry::nAmp];
  Int_t binNum;
  if(BS->dimensions==1) {
    binNum = BS->HashBinNum(0);

    kindepGr = kindepMap.at(binNum);
    kindepCanv->cd();
    DrawKinDepGraph(kindepGr,BS,0);

    multiGr = multiMap.at(binNum);
    if(includeOAonMultiGr) multiGr->Add(kindepGr);
    for(int aa=0; aa<N_AMP; aa++) {
      RFkindepGr[aa] = RFkindepMap[aa].at(binNum);
      RFkindepCanv[aa]->cd();
      DrawKinDepGraph(RFkindepGr[aa],BS,0);
      RFkindepGrClone[aa] = ShiftGraph(RFkindepGr[aa],aa+1);
      multiGr->Add(RFkindepGrClone[aa]);
    };

    chindfGr = chindfMap.at(binNum);
    chindfCanv->cd();
    DrawSimpleGraph(chindfGr,BS,0);

    rellumGr = rellumMap.at(binNum);
    rellumCanv->cd();
    DrawSimpleGraph(rellumGr,BS,0,false);

    for(int b0=0; b0<NB[0]; b0++) {
      binNum = BS->HashBinNum(b0);
      A = asymMap.at(binNum);
      asymModCanv->cd(b0+1);
      if(!(A->oa2d)) DrawAsymGr(A->asymGr);
      else DrawAsymGr2(A->asymGr2);
      modDistCanv->cd(b0+1);
      if(!(A->oa2d)) A->modDist->Draw();
      else A->modDist2->Draw("colz");
      if(A->oa2d) {
        asymModHist2Canv->cd(b0+1);
        A->asymGr2hist->Draw("colz");
      };
    };
  }
  else if(BS->dimensions==2) {
    for(int b1=0; b1<NB[1]; b1++) {
      pad = b1+1;
      binNum = BS->HashBinNum(0,b1);

      kindepGr = kindepMap.at(binNum);
      kindepCanv->cd(pad);
      DrawKinDepGraph(kindepGr,BS,0);

      multiGr = multiMap.at(binNum);
      if(includeOAonMultiGr) multiGr->Add(kindepGr);

      for(int aa=0; aa<N_AMP; aa++) {
        RFkindepGr[aa] = RFkindepMap[aa].at(binNum);
        RFkindepCanv[aa]->cd(pad);
        DrawKinDepGraph(RFkindepGr[aa],BS,0);
        RFkindepGrClone[aa] = ShiftGraph(RFkindepGr[aa],aa+1);
        multiGr->Add(RFkindepGrClone[aa]);
      };

      chindfGr = chindfMap.at(binNum);
      chindfCanv->cd(pad);
      DrawSimpleGraph(chindfGr,BS,0);

      rellumGr = rellumMap.at(binNum);
      rellumCanv->cd(pad);
      DrawSimpleGraph(rellumGr,BS,0,false);

      for(int b0=0; b0<NB[0]; b0++) {
        binNum = BS->HashBinNum(b0,b1);
        A = asymMap.at(binNum);
        asymModCanv->cd(b0*NB[1]+b1+1);
        if(!(A->oa2d)) DrawAsymGr(A->asymGr);
        else DrawAsymGr2(A->asymGr2);
        modDistCanv->cd(b0*NB[1]+b1+1);
        if(!(A->oa2d)) A->modDist->Draw();
        else A->modDist2->Draw();
        if(A->oa2d) {
          asymModHist2Canv->cd(b0*NB[1]+b1+1);
          A->asymGr2hist->Draw("colz");
        };
      };
    };
  }
  else if(BS->dimensions==3) {
    for(int b1=0; b1<NB[1]; b1++) {
      for(int b2=0; b2<NB[2]; b2++) {
        pad = b1*NB[2]+b2+1;
        binNum = BS->HashBinNum(0,b1,b2);

        kindepGr = kindepMap.at(binNum);
        kindepCanv->cd(pad);
        DrawKinDepGraph(kindepGr,BS,0);

        multiGr = multiMap.at(binNum);
        if(includeOAonMultiGr) multiGr->Add(kindepGr);

        for(int aa=0; aa<N_AMP; aa++) {
          RFkindepGr[aa] = RFkindepMap[aa].at(binNum);
          RFkindepCanv[aa]->cd(pad);
          DrawKinDepGraph(RFkindepGr[aa],BS,0);
          RFkindepGrClone[aa] = ShiftGraph(RFkindepGr[aa],aa+1);
          multiGr->Add(RFkindepGrClone[aa]);
        };

        chindfGr = chindfMap.at(binNum);
        chindfCanv->cd(pad);
        DrawSimpleGraph(chindfGr,BS,0);

        rellumGr = rellumMap.at(binNum);
        rellumCanv->cd(pad);
        DrawSimpleGraph(rellumGr,BS,0,false);
      };
    };
  };


  // build multiGr canvases
  TCanvas * multiGrCanv;
  TLegend * multiLeg;
  TString legText;
  TString multiGrCanvN;
  TObjArray * multiGrCanvArr;
  Int_t nRows = ( N_AMP + 1 - (includeOAonMultiGr?0:1) )/4 + 1;
  for(Int_t bn : BS->binVec) {
    A = asymMap.at(bn);
    if(A->B[0] == 0) {

      multiGr = multiMap.at(bn);
      multiGrCanvN = multiGr->GetName();
      multiGrCanvN.ReplaceAll("multiGr","multiGrCanv");
      multiGrCanv = new TCanvas(multiGrCanvN,multiGrCanvN,1600,nRows*400);
      multiGrCanv->Divide(4,nRows);

      multiLeg = new TLegend(0.1,0.1,0.9,0.9);

      modu->enablePW = A->fitPW;

      for(int aa=0; aa<N_AMP; aa++) {
        multiGrCanv->cd(aa+1);
        multiGrCanv->GetPad(aa+1)->SetGrid(0,1);
        RFkindepGr[aa] = RFkindepMap[aa].at(bn);
        RFkindepGr[aa]->Draw("LAPE");
        legText = modu->StateTitle(A->fitTw[aa],A->fitL[aa],A->fitM[aa]);
        legText += ": ";
        legText += modu->ModulationTitle(A->fitTw[aa],A->fitL[aa],A->fitM[aa]);
        multiLeg->AddEntry(RFkindepGr[aa],legText,"PLE");
      };
      
      if(includeOAonMultiGr) {
        kindepGr = kindepMap.at(bn);
        multiGrCanv->cd(N_AMP+1);
        multiGrCanv->GetPad(N_AMP+1)->SetGrid(0,1);
        kindepGr->Draw("LAPE");
        legText = A->oaModulationTitle;
        legText += " one-amp result";
        multiLeg->AddEntry(kindepGr,legText,"PLE");
      };

      multiGrCanv->cd(4*nRows);
      multiLeg->Draw();

      multiGrCanvArr = new TObjArray();
      multiGrCanvArr->AddLast(multiGrCanv);

    };
  };


  // sum distributions (for showing bin boundaries)
  TH1D * ivFullDist1;
  TH2D * ivFullDist2;
  TH3D * ivFullDist3;
  TH1D * modFullDist;
  TH2D * modFullDist2; // for 2d modulation
  TH2D * IVvsModFullDist;

  if(BS->dimensions==1) {
    for(int b0=0; b0<NB[0]; b0++) {
      binNum = BS->HashBinNum(b0);
      A = asymMap.at(binNum);
      if(b0==0) {
        ivFullDist1 = (TH1D*)(A->ivDist1)->Clone();
        SetCloneName(ivFullDist1);
        if(!(A->oa2d)) {
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
        if(!(A->oa2d)) {
          IVvsModFullDist->Add(A->IVvsModDist);
          modFullDist->Add(A->modDist);
        } else {
          modFullDist2->Add(A->modDist2);
        };
      };
    };
  }
  else if(BS->dimensions==2) {
    for(int b1=0; b1<NB[1]; b1++) {
      for(int b0=0; b0<NB[0]; b0++) {
        binNum = BS->HashBinNum(b0,b1);
        A = asymMap.at(binNum);
        if(b0==0 && b1==0) {
          ivFullDist2 = (TH2D*)(A->ivDist2)->Clone();
          SetCloneName(ivFullDist2);
          if(!(A->oa2d)) {
            modFullDist = (TH1D*)(A->modDist)->Clone();
            SetCloneName(modFullDist);
          } else {
            modFullDist2 = (TH2D*)(A->modDist2)->Clone();
            SetCloneName(modFullDist2);
          };
        } else {
          ivFullDist2->Add(A->ivDist2);
          if(!(A->oa2d)) modFullDist->Add(A->modDist);
          else modFullDist2->Add(A->modDist2);
        };
      };
    };
  }
  else if(BS->dimensions==3) {
    for(int b2=0; b2<NB[2]; b2++) {
      for(int b1=0; b1<NB[1]; b1++) {
        for(int b0=0; b0<NB[0]; b0++) {
          binNum = BS->HashBinNum(b0,b1,b2);
          A = asymMap.at(binNum);
          if(b0==0 && b1==0 && b2==0) {
            ivFullDist3 = (TH3D*)(A->ivDist3)->Clone();
            SetCloneName(ivFullDist3);
            if(!(A->oa2d)) {
              modFullDist = (TH1D*)(A->modDist)->Clone();
              SetCloneName(modFullDist);
            } else {
              modFullDist2 = (TH2D*)(A->modDist2)->Clone();
              SetCloneName(modFullDist2);
            };
          } else {
            ivFullDist3->Add(A->ivDist3);
            if(!(A->oa2d)) modFullDist->Add(A->modDist);
            else modFullDist2->Add(A->modDist2);
          };
        };
      };
    };
  };


  // write output to asymFile
  asymFile->cd();
  // -- Asymmetry objects
  printf("--- write Asymmetry objects\n");
  for(Int_t bn : BS->binVec) {
    A = asymMap.at(bn);
    A->StreamData(asymFile);
  };

  // -- "full" distributions
  if(BS->dimensions==1) {
    ivFullDist1->Write();
    if(!(A->oa2d)) IVvsModFullDist->Write();
  } else if(BS->dimensions==2) {
    ivFullDist2->Write();
  } else if(BS->dimensions==3) {
    ivFullDist3->Write();
  };
  if(!(A->oa2d)) modFullDist->Write();
  else modFullDist2->Write();

  // -- asymmetries and kindep graphs
  printf("--- write kinematic-dependent asymmetries\n");
  for(Int_t bn : BS->binVec) {
    A = asymMap.at(bn);
    A->PrintSettings();

    // first write out the asymmetry vs. modulation graphs
    if(!(A->oa2d)) A->asymGr->Write();
    else A->asymGr2->Write();

    // then write out the kindepGr after writing out all the
    // relevant asymmetry vs. modulation graphs
    if(A->B[0] + 1 == NB[0]) {
      binNum = BS->HashBinNum(A->B[0], A->B[1], A->B[2]);
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

  // write multiGr canvases
  multiGrCanvArr->Write("multiGrCanvArr",TObject::kSingleKey);


  kindepCanv->Write();
  for(int aa=0; aa<N_AMP; aa++) RFkindepCanv[aa]->Write();
  chindfCanv->Write();
  chisqDist->Write();
  rellumCanv->Write();
  if(BS->dimensions==1 || BS->dimensions==2) {
    asymModCanv->Write();
    if(A->oa2d) asymModHist2Canv->Write();
    modDistCanv->Write();
  };


  // -- RooFit results
  TCanvas * rfCanv[Asymmetry::nAmp];
  TString rfCanvName[Asymmetry::nAmp];

  for(Int_t bn : BS->binVec) {
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


  asymFile->Close();
  catFile->Close();
  printf("--- end %s\n",argv[0]);
  return 0;

};



void DrawKinDepGraph(TGraphErrors * g_, Binning * B_, Int_t d_) {
  Int_t v_ = B_->ivVar[d_];

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



void DrawSimpleGraph(TGraphErrors * g_, Binning * B_, Int_t d_, Bool_t setRange) {
  Int_t v_ = B_->ivVar[d_];

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
      retGr->SetLineColor(N_AMP==1?kGray+3:kGray+1); 
      //retGr->SetLineStyle(2);
      retGr->SetMarkerStyle(kFullCircle);
      break;
    case 2:
      retGr->SetLineColor(kRed); 
      //retGr->SetLineStyle(1);
      retGr->SetMarkerStyle(kFullCircle);
      break;
    case 3:
      retGr->SetLineColor(kAzure+1);
      //retGr->SetLineStyle(3);
      retGr->SetMarkerStyle(kFullCircle);
      break;
    case 4:
      retGr->SetLineColor(kViolet+1);
      retGr->SetMarkerStyle(kFullCircle);
      break;
    case 5:
      retGr->SetLineColor(kGreen+1);
      retGr->SetMarkerStyle(kFullCircle);
      break;
    case 6:
      retGr->SetLineColor(kCyan);
      retGr->SetMarkerStyle(kFullCircle);
      break;
    case 7:
      retGr->SetLineColor(kOrange);
      retGr->SetMarkerStyle(kFullCircle);
      break;
    case 8:
      retGr->SetLineColor(kOrange-7);
      retGr->SetMarkerStyle(kFullCircle);
      break;
    default: retGr->SetLineColor(kGray);
  };


  retGr->SetMarkerColor(kBlack);
  retGr->SetLineWidth(2);
  retGr->SetMarkerSize(1.3);

  gr->SetLineColor(retGr->GetLineColor());
  gr->SetLineStyle(retGr->GetLineStyle());
  gr->SetMarkerStyle(retGr->GetMarkerStyle());
  gr->SetMarkerColor(kBlack);
  gr->SetLineWidth(2);
  gr->SetMarkerSize(1.3);

  return retGr;
};
