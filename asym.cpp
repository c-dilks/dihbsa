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

// DihBsa
#include "Constants.h"
#include "DIS.h"
#include "Trajectory.h"
#include "Dihadron.h"
#include "EventTree.h"
#include "Binning.h"
#include "Asymmetry.h"

Int_t GetBinNum(Int_t bin0, Int_t bin1=-1, Int_t bin2=-1);
void DrawKinDepGraph(TGraphErrors * g_, Binning * B_, Int_t v_);
void DrawSimpleGraph(TGraphErrors * g_, Binning * B_, Int_t v, Bool_t setRange_=true);
void DrawAsymGr(TGraphErrors * g_);
void DrawAsymGr2(TGraph2DErrors * g_);
void SetCloneName(TH1 * clone_);

Int_t pairType;
TString inDir;
Int_t whichModulation;
Int_t dimensions;
Int_t ivType;
Int_t whichPhiR;
Bool_t batchMode;

TString dihTitle;

int main(int argc, char** argv) {
   
   gStyle->SetOptFit(1);


   // ARGUMENTS
   pairType = EncodePairType(kPip,kPim);
   inDir = "outroot";
   whichModulation = 0; // see src/Asymmetry.h
   dimensions = 1; // number of dimensions 
   ivType = 1; // which variables to bin in (see below)
   whichPhiR = 3; // 1:phiRq  2:phiRp_r  3:phiRp // Alessandro prefers 3:phiRp
   batchMode = 0; // if true, renames output files and prints pngs

   if(argc>1) inDir =           TString(argv[1]);
   if(argc>2) pairType =        (Int_t) strtof(argv[2],NULL);
   if(argc>3) whichModulation = (Int_t) strtof(argv[3],NULL);
   if(argc>4) dimensions =      (Int_t) strtof(argv[4],NULL);
   if(argc>5) ivType =          (Int_t) strtof(argv[5],NULL);
   if(argc>6) whichPhiR =       (Int_t) strtof(argv[6],NULL);
   if(argc>7) batchMode =       (Bool_t) strtof(argv[7],NULL);

   printf("pairType = 0x%x\n",pairType);
   printf("inDir = %s\n",inDir.Data());
   printf("whichModulation = %d\n",whichModulation);
   printf("dimensions = %d\n",dimensions);
   printf("ivType = %d\n",ivType);
   printf("whichPhiR = %d\n",whichPhiR);
   printf("batchMode = %d\n",(Int_t)batchMode);
   printf("\n");

   Binning * BS = new Binning(pairType); // instantiate binning scheme 
   Asymmetry * A; // Asymmetry pointer
   Int_t whichHad[2];
   DecodePairType(pairType,whichHad[qA],whichHad[qB]);
   dihTitle = PairTitle(whichHad[qA],whichHad[qB]);

   // help printout
   if(argc==1) {
     fprintf(stderr,
       "\nUSAGE: %s [inDir] [pairType] [whichModulation] [dimensions] [ivType] [whichPhiR] [batchMode]\n",
       argv[0]);
     printf("\n- inDir: directory of ROOT files to analyse\n");
     printf("\n- pairType: hadron pair type (run PrintEnumerators.C)\n");
     printf("\n- whichModulation:\n");
     for(int m=0; m<Asymmetry::nMod; m++) {
       A = new Asymmetry(BS,m,-10000);
       printf("   %d = %s =  %s\n",m,
         (A->ModulationName).Data(),(A->ModulationTitle).Data());
     };
     printf("\n- dimensions: for multi-dimensional asymmetry analysis in 1, 2, or 3-D\n");
     printf("\n- ivType: 3-digit number, one for each dimension, ");
     printf("where the digits represent:\n");
     for(int i=0; i<Binning::nIV; i++) {
       printf("   %d = %s\n",i,(BS->IVtitle[i]).Data());
     };
     printf("\n- whichPhiR:\n");
     printf("   1 = PhiRq: from R_perp via vector rejection\n");
     printf("   2 = PhiRp_r: from R_T via vector rejection\n");
     printf("   3 = PhiRp: from R_T via k_T relation (PREFERRED)\n");
     printf("\n- batchMode: boolean, where if true, renames output root file\n");
     printf("  and prints pngs\n");
     printf("\n");
     return 0;
   };



   // ivType:
   // -- if dimensions==1, all asymetries will be plotted against one independent
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
       ivVar[0] = ivType;
       break;
     case 2:
       ivVar[0] = ivType / 10; 
       ivVar[1] = ivType % 10;
       break;
     case 3:
       ivVar[0] = ivType / 100;
       ivVar[1] = ( ivType / 10 ) % 10;
       ivVar[2] = ivType % 10;
       break;
     default:
       fprintf(stderr,"ERROR: bad number of dimensions\n");
       return 0;
   };


   // check IV enumerators and get number of bins for each IV
   Int_t NB[3];
   for(int d=0; d<dimensions; d++) {
     printf("ivVar[%d] = %d\n",d,ivVar[d]);
     if(!(BS->ValidIV(ivVar[d]))) {
       fprintf(stderr,"ERROR: this IV is unknown\n");
       return 0;
     };
     NB[d] = BS->nBins[ivVar[d]];
   };


   // print which IV will be analyzed
   A = new Asymmetry(BS,whichModulation,-10000);
   TString modN = A->ModulationName;
   printf("--------> Analysing %s asymmetries vs. %s ",
     modN.Data(),(BS->IVname[ivVar[0]]).Data());
   if(dimensions>=2)
     printf("in bins of %s ",(BS->IVname[ivVar[1]]).Data());
   if(dimensions>=3)
     printf("and %s ",(BS->IVname[ivVar[2]]).Data());
   printf("\n\n");
   

   // set output file
   TString outfileName = "spin";
   if(batchMode) {
     for(int d=0; d<dimensions; d++) 
       outfileName = outfileName + "_" + BS->IVname[ivVar[d]];
     outfileName = outfileName + "_" + modN;
   };
   outfileName = outfileName + ".root";
   printf("outfileName = %s\n",outfileName.Data());
   TFile * outfile = new TFile(outfileName,"RECREATE");


   EventTree * ev = new EventTree(TString(inDir+"/*.root"),pairType);


   
   // instantiate Asymmetry pointers 
   std::vector<Asymmetry*> asymVec; // vector of Asymmetry instances, one for each bin
                                    // (iterator invalidation is not protected
                                    // against, but access may be faster than a map)
   std::map<Int_t, Asymmetry*> asymMap; // map of asymVec entry number (iterator) 
                                        // and Asymmetry instance (this may be slower
                                        // than asymVec, but protects against iterator
                                        // invalidation)
   std::map<int,int> binMap; // 3-digit bin number -> asym iterator
                                 // -- use GetBinNum(...) to convert a set of bin
                                 //    numbers to the 3-digit bin number; then use this
                                 //    map to get the index of the corresponding
                                 //    Asymmetry pointer in asymVec
   std::map<int,TGraphErrors*> kindepMap; // 3-digit bin number -> kinematic-dependent
                                          // -- use GetBinNum(...) to convert a set of
                                          // bin numbers to the 3-digit bin number; then
                                          // use this map to get the index of the
                                          // corresponding kinematic-dependent asymmetry
                                          // graph
   std::map<int,TGraphErrors*> chindfMap; // kindepMap for chindfGr
   std::map<int,TGraphErrors*> rellumMap; // kindepMap for rellumGr

   TGraphErrors * kindepGr;
   TGraphErrors * RFkindepGr;
   TGraphErrors * chindfGr;
   TGraphErrors * rellumGr;
   TString grTitle,grName;

   Int_t binNum;
   Int_t bcnt = 0;
   if(dimensions == 1) {

     for(int b=0; b<NB[0]; b++) {

       A = new Asymmetry(BS, whichModulation, 1, ivVar[0], b);
       if(!(A->success)) return 0;
       asymVec.push_back(A);
       binNum = GetBinNum(b);
       binMap.insert(std::pair<int,int>(binNum,bcnt));
       asymMap.insert(std::pair<Int_t, Asymmetry*>(binNum,A));
       bcnt++;
       if(b==0) {
         grTitle = Form("%s %s asymmetry vs. %s;%s",
           dihTitle.Data(),(A->ModulationTitle).Data(),
           (BS->IVtitle[ivVar[0]]).Data(),(BS->IVtitle[ivVar[0]]).Data());
         grName = Form("kindep_%s",(BS->IVname[ivVar[0]]).Data());
         kindepGr = new TGraphErrors();
         kindepGr->SetName(grName);
         kindepGr->SetTitle(grTitle);

         RFkindepGr = new TGraphErrors();
         RFkindepGr->SetName(TString("RF"+grName));
         RFkindepGr->SetTitle(grTitle);
         // aqui

         grTitle = "#chi^{2}/NDF of " + grTitle;
         grName.ReplaceAll("kindep","chindf");
         chindfGr = new TGraphErrors();
         chindfGr->SetName(grName);
         chindfGr->SetTitle(grTitle);

         grTitle = Form("relative luminosity vs. %s;%s",
           (BS->IVtitle[ivVar[0]]).Data(),(BS->IVtitle[ivVar[0]]).Data());
         grName.ReplaceAll("chindf","rellum");
         rellumGr = new TGraphErrors();
         rellumGr->SetName(grName);
         rellumGr->SetTitle(grTitle);
       };

       kindepMap.insert(std::pair<Int_t,TGraphErrors*>(binNum,kindepGr));
       chindfMap.insert(std::pair<Int_t,TGraphErrors*>(binNum,chindfGr));
       rellumMap.insert(std::pair<Int_t,TGraphErrors*>(binNum,rellumGr));

     };
   }
   else if(dimensions == 2) {
     for(int b1=0; b1<NB[1]; b1++) {
       for(int b0=0; b0<NB[0]; b0++) {

         A = new Asymmetry(BS, whichModulation, 2, 
           ivVar[0], b0,
           ivVar[1], b1
         );
         if(!(A->success)) return 0;
         asymVec.push_back(A);
         binNum = GetBinNum(b0,b1);
         binMap.insert(std::pair<int,int>(binNum,bcnt));
         asymMap.insert(std::pair<Int_t, Asymmetry*>(binNum,A));
         bcnt++;

         if(b0==0) {
           grTitle = Form("%s %s asymmetry vs. %s :: %s;%s",
             dihTitle.Data(),
             (A->ModulationTitle).Data(),(BS->IVtitle[ivVar[0]]).Data(),
             (BS->GetBoundStr(ivVar[1],b1)).Data(),
             (BS->IVtitle[ivVar[0]]).Data());
           grName = Form("kindep_%s_bin_%s%d",(BS->IVname[ivVar[0]]).Data(),
             (BS->IVname[ivVar[1]]).Data(),b1);
           kindepGr = new TGraphErrors();
           kindepGr->SetName(grName);
           kindepGr->SetTitle(grTitle);

           grTitle = "#chi^{2}/NDF of " + grTitle;
           grName.ReplaceAll("kindep","chindf");
           chindfGr = new TGraphErrors();
           chindfGr->SetName(grName);
           chindfGr->SetTitle(grTitle);

           grTitle = Form("relative luminosity vs. %s :: %s;%s",
             (BS->IVtitle[ivVar[0]]).Data(),
             (BS->GetBoundStr(ivVar[1],b1)).Data(),(BS->IVtitle[ivVar[0]]).Data());
           grName.ReplaceAll("chindf","rellum");
           rellumGr = new TGraphErrors();
           rellumGr->SetName(grName);
           rellumGr->SetTitle(grTitle);
         };

         kindepMap.insert(std::pair<Int_t,TGraphErrors*>(binNum,kindepGr));
         chindfMap.insert(std::pair<Int_t,TGraphErrors*>(binNum,chindfGr));
         rellumMap.insert(std::pair<Int_t,TGraphErrors*>(binNum,rellumGr));

       };
     };
   }
   else if(dimensions == 3) {
     for(int b2=0; b2<NB[2]; b2++) {
       for(int b1=0; b1<NB[1]; b1++) {
         for(int b0=0; b0<NB[0]; b0++) {

           A = new Asymmetry(BS, whichModulation, 3, 
             ivVar[0], b0,
             ivVar[1], b1,
             ivVar[2], b2
           );
           if(!(A->success)) return 0;
           asymVec.push_back(A);
           binNum = GetBinNum(b0,b1,b2);
           binMap.insert(std::pair<int,int>(binNum,bcnt));
           asymMap.insert(std::pair<Int_t, Asymmetry*>(binNum,A));
           bcnt++;

           if(b0==0) {
             grTitle = Form("%s %s asymmetry vs. %s :: %s, %s;%s",
               dihTitle.Data(),
               (A->ModulationTitle).Data(),(BS->IVtitle[ivVar[0]]).Data(),
               (BS->GetBoundStr(ivVar[1],b1)).Data(),
               (BS->GetBoundStr(ivVar[2],b2)).Data(),
               (BS->IVtitle[ivVar[0]]).Data());
             grName = Form("kindep_%s_bin_%s%d_%s%d",(BS->IVname[ivVar[0]]).Data(),
               (BS->IVname[ivVar[1]]).Data(),b1,(BS->IVname[ivVar[2]]).Data(),b2);
             kindepGr = new TGraphErrors();
             kindepGr->SetName(grName);
             kindepGr->SetTitle(grTitle);

             grTitle = "#chi^{2}/NDF of " + grTitle;
             grName.ReplaceAll("kindep","chindf");
             chindfGr = new TGraphErrors();
             chindfGr->SetName(grName);
             chindfGr->SetTitle(grTitle);

             grTitle = Form("relative luminosity vs. %s :: %s, %s;%s",
               (BS->IVtitle[ivVar[0]]).Data(),
               (BS->GetBoundStr(ivVar[1],b1)).Data(),
               (BS->GetBoundStr(ivVar[2],b2)).Data(),
               (BS->IVtitle[ivVar[0]]).Data());
             grName.ReplaceAll("chindf","rellum");
             rellumGr = new TGraphErrors();
             rellumGr->SetName(grName);
             rellumGr->SetTitle(grTitle);
           };

           kindepMap.insert(std::pair<Int_t,TGraphErrors*>(binNum,kindepGr));
           chindfMap.insert(std::pair<Int_t,TGraphErrors*>(binNum,chindfGr));
           rellumMap.insert(std::pair<Int_t,TGraphErrors*>(binNum,rellumGr));

         };
       };
     };
   };


   // overall summary plots
   TH1F * chisqDist = new TH1F("chisqDist","#chi^{2} distribution",100,0,20);


   // EVENT LOOP -------------------------------------------
   printf("begin loop through %lld events...\n",ev->ENT);
   Bool_t filled;
   for(int i=0; i<ev->ENT; i++) {

     ev->GetEvent(i);

     if(ev->Valid()) {
       
       // fill asymmetry plots; Asymmetry::FillPlots() checks the bin,
       // and fills plots if it's the correct bin
       for(std::vector<Asymmetry*>::iterator it = asymVec.begin(); 
         it!=asymVec.end(); ++it
       ) {
         
         A = *it;

         switch(whichPhiR) {
           case 1:
             A->PhiR = ev->PhiRq;
             break;
           case 2:
             A->PhiR = ev->PhiRp_r;
             break;
           case 3:
             A->PhiR = ev->PhiRp;
             break;
           default:
             fprintf(stderr,"ERROR: invalid whichPhiR\n");
             return 0;
         };

         A->Mh = ev->Mh;
         A->x = ev->x;
         A->z = ev->Zpair;
         A->eSpin = ev->helicity;
         A->pSpin = 0;
         A->PhiH = ev->PhiH;
         A->PhPerp = ev->PhPerp;
         A->theta = ev->theta;

         filled = A->FillPlots(); 

         //if(filled && A->debug) ev->PrintEvent();
       };

     };
   };
   // end event loop -------------------------------------------



   // compute asymmetries
   Float_t asymValue,asymError;
   Float_t kinValue,kinError;
   Float_t chisq,ndf;
   printf("--- calculate asymmetries\n");
   for(std::vector<Asymmetry*>::iterator it = asymVec.begin(); 
     it!=asymVec.end(); ++it
   ) {

     A = *it;
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
       binNum = GetBinNum(A->B[0], A->B[1], A->B[2]);

       kindepGr = kindepMap.at(binNum);
       kindepGr->SetPoint(A->B[0],kinValue,asymValue);
       kindepGr->SetPointError(A->B[0],kinError,asymError);

       chindfGr = chindfMap.at(binNum);
       chindfGr->SetPoint(A->B[0],kinValue,chisq/ndf);

       rellumGr = rellumMap.at(binNum);
       rellumGr->SetPoint(A->B[0],kinValue,A->rellum);
       rellumGr->SetPointError(A->B[0],0,A->rellumErr);

     };

   };


   //gStyle->SetOptFit(1); // (better to put this in your ~/.rootlogon.C file)

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




   // -- add objects to canvases
   Int_t pad;
   if(dimensions==1) {
     binNum = GetBinNum(0);
     kindepGr = kindepMap.at(binNum);
     kindepCanv->cd();
     DrawKinDepGraph(kindepGr,BS,ivVar[0]);

     chindfGr = chindfMap.at(binNum);
     chindfCanv->cd();
     DrawSimpleGraph(chindfGr,BS,ivVar[0]);

     rellumGr = rellumMap.at(binNum);
     rellumCanv->cd();
     DrawSimpleGraph(rellumGr,BS,ivVar[0],false);

     for(int b0=0; b0<NB[0]; b0++) {
       binNum = GetBinNum(b0);
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
       binNum = GetBinNum(0,b1);
       kindepGr = kindepMap.at(binNum);
       kindepCanv->cd(pad);
       DrawKinDepGraph(kindepGr,BS,ivVar[0]);
       
       chindfGr = chindfMap.at(binNum);
       chindfCanv->cd(pad);
       DrawSimpleGraph(chindfGr,BS,ivVar[0]);
       
       rellumGr = rellumMap.at(binNum);
       rellumCanv->cd(pad);
       DrawSimpleGraph(rellumGr,BS,ivVar[0],false);

       for(int b0=0; b0<NB[0]; b0++) {
         binNum = GetBinNum(b0,b1);
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
         binNum = GetBinNum(0,b1,b2);
         kindepGr = kindepMap.at(binNum);
         kindepCanv->cd(pad);
         DrawKinDepGraph(kindepGr,BS,ivVar[0]);

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
       binNum = GetBinNum(b0);
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
         binNum = GetBinNum(b0,b1);
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
           binNum = GetBinNum(b0,b1,b2);
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
   // -- Asymmetry objects
   printf("--- write Asymmetry objects\n");
   for(std::vector<Asymmetry*>::iterator it = asymVec.begin(); 
     it!=asymVec.end(); ++it
   ) {

     A = *it;

     printf("writing plots for:\n");
     A->PrintSettings();

     switch(dimensions) {
       case 1: A->ivDist1->Write(); break;
       case 2: A->ivDist2->Write(); break;
       case 3: A->ivDist3->Write(); break;
     };

     if(!(A->asym2d)) {
       A->modDist->Write();
       for(Int_t m=0; m<Asymmetry::nModBins; m++) A->modBinDist[m]->Write();
       if(dimensions==1) A->IVvsModDist->Write();
       for(Int_t s=0; s<nSpin; s++) A->aziDist[s]->Write();
     } else {
       A->modDist2->Write();
       for(Int_t mmH=0; mmH<Asymmetry::nModBins2; mmH++) {
         for(Int_t mmR=0; mmR<Asymmetry::nModBins2; mmR++) {
           A->modBinDist2[mmH][mmR]->Write();
         };
       };
       for(Int_t s=0; s<nSpin; s++) A->aziDist2[s]->Write();
     };
   };

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
   for(std::vector<Asymmetry*>::iterator it = asymVec.begin(); 
     it!=asymVec.end(); ++it
   ) {
     A = *it;
     A->PrintSettings();

     // first write out the asymmetry vs. modulation graphs
     if(!(A->asym2d)) A->asymGr->Write();
     else A->asymGr2->Write();

     // then write out the kindepGr *after* writing out all the
     // relevant asymmetry vs. modulation graphs
     if(A->B[0] + 1 == NB[0]) {
       binNum = GetBinNum(A->B[0], A->B[1], A->B[2]);
       kindepGr = kindepMap.at(binNum);
       kindepGr->Write();
     };
   };

   kindepCanv->Write();
   chindfCanv->Write();
   chisqDist->Write();
   rellumCanv->Write();
   if(dimensions==1 || dimensions==2) {
     asymModCanv->Write();
     if(A->asym2d) asymModHist2Canv->Write();
     modDistCanv->Write();
   };


   // TCanvas for RooFit result
   TCanvas * rfCanv[nSpin];
   TGraphErrors * rfAsymGr = new TGraphErrors();
   rfAsymGr->SetName("roofitResult");
   rfAsymGr->SetTitle("roofitResult");
   rfAsymGr->SetMarkerStyle(kFullCircle);
   rfAsymGr->SetLineColor(kBlue);
   Int_t rfAsymGrCnt = 0;
   TString rfCanvName[nSpin];
   for(std::vector<Asymmetry*>::iterator it = asymVec.begin(); 
     it!=asymVec.end(); ++it
   ) {
     A = *it;
     if(A->roofitter) {
       
       for(int ss=0; ss<nSpin; ss++) {
         rfCanvName[ss] = "RF_spin" + SpinName(ss) + "_" + A->asymGr->GetName();
         rfCanv[ss] = new TCanvas(rfCanvName[ss],rfCanvName[ss],800,800);
         A->rfPhiRplot[ss]->Draw();
         rfCanv[ss]->Write();
       };

       A->PrintSettings();
       printf("rfAR = %f +- %f\n",A->rfAR->getVal(),A->rfAR->getError());

       rfAsymGr->SetPoint(rfAsymGrCnt,rfAsymGrCnt+1,A->rfAR->getVal());
       rfAsymGr->SetPointError(rfAsymGrCnt,0,A->rfAR->getError());
       rfAsymGrCnt++;

     };
   };
   rfAsymGr->Write();




   // print modDist boundaries (used for determining modDist boundaries
   // which depend on kinematics, for g1perp modulation PhPerp/Mh-scaling test)
   Float_t mbound;
   if(dimensions==1 && whichModulation==Asymmetry::scaleSinPhiHR) {
     printf("MBOUND\n");
     printf("if(v_==v%s) {\n",(BS->IVname[ivVar[0]]).Data());
     for(int b=0; b<NB[0]; b++) {
       binNum = GetBinNum(b);
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



   // print images
   TString pngName;
   if(batchMode) {
     pngName = Form("%s.png",kindepCanv->GetName());
     kindepCanv->Print(pngName,"png");
     pngName = Form("%s.png",chindfCanv->GetName());
     chindfCanv->Print(pngName,"png");
     pngName = Form("%s.png",rellumCanv->GetName());
     rellumCanv->Print(pngName,"png");
     if(dimensions==1 || dimensions==2) {
       pngName = Form("%s.png",asymModCanv->GetName());
       asymModCanv->Print(pngName,"png"); 
       pngName = Form("%s.png",modDistCanv->GetName());
       modDistCanv->Print(pngName,"png"); 
     };
   };

   outfile->Close();

   printf("--- end %s\n",argv[0]);

};


Int_t GetBinNum(Int_t bin0, Int_t bin1, Int_t bin2) {
  Int_t retval = bin0;
  if(bin1>=0) retval += 10 * bin1;
  if(bin2>=0) retval += 100 * bin2;
  return retval;
};


void DrawKinDepGraph(TGraphErrors * g_, Binning * B_, Int_t v_) {

  g_->Draw("APE"); // draw once, so we can then format it

  g_->SetLineColor(B_->GetColor(v_));
  g_->SetLineWidth(2);

  g_->SetMarkerStyle(kFullCircle);
  g_->SetMarkerColor(kBlack);
  g_->SetMarkerSize(1.3);

  // set vertical axis range (it is overridden if the plot's vertical range
  // is larger than the desired range)
  Float_t yMin = -0.07;
  Float_t yMax = 0.1;
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

  //g_->SetLineColor(B_->GetColor(v_));
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
