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


int main(int argc, char** argv) {

   // ARGUMENTS
   TString inDir = "outroot";
   Int_t whichModulation = 0; // see src/Asymmetry.h
   Int_t dimensions = 1; // number of dimensions 
   Int_t ivType = 0; // which variables to bin in (see below)
   Int_t whichPhiR = 3; // 1:phiRq  2:phiRp_r  3:phiRp // Alessandro prefers 3:phiRp

   if(argc>1) inDir = TString(argv[1]);
   if(argc>2) whichModulation = (Int_t)strtof(argv[2],NULL);
   if(argc>3) dimensions = (Int_t)strtof(argv[3],NULL);
   if(argc>3) ivType = (Int_t)strtof(argv[4],NULL);
   if(argc>5) whichPhiR = (Int_t)strtof(argv[5],NULL);

   printf("inDir = %s\n",inDir.Data());
   printf("whichModulation = %d\n",whichModulation);
   printf("dimensions = %d\n",dimensions);
   printf("ivType = %d\n",ivType);
   printf("whichPhiR = %d\n",whichPhiR);
   printf("\n");

   // instantiate binning scheme
   Binning * BS = new Binning();

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
   printf("--------> Analysing asymmetries vs. %s ",(BS->IVname[ivVar[0]]).Data());
   if(dimensions>=2)
     printf("in bins of %s ",(BS->IVname[ivVar[1]]).Data());
   if(dimensions>=3)
     printf("and %s ",(BS->IVname[ivVar[2]]).Data());
   printf("\n\n");
   

   // set output file
   TFile * outfile = new TFile("spin.root","RECREATE");
   EventTree * ev = new EventTree(TString(inDir+"/*.root"));

   
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

   Asymmetry * A;
   TGraphErrors * kindepGr;
   TString grTitle,grName;

   Int_t binNum;
   Int_t bcnt = 0;
   if(dimensions == 1) {

     for(int b=0; b<NB[0]; b++) {

       A = new Asymmetry(BS, whichModulation, 1, ivVar[0], b);
       asymVec.push_back(A);
       binNum = GetBinNum(b);
       binMap.insert(std::pair<int,int>(binNum,bcnt));
       asymMap.insert(std::pair<Int_t, Asymmetry*>(binNum,A));
       bcnt++;

       if(b==0) {
         grTitle = Form("%s asymmetry vs. %s",
           A->ModulationTitle.Data(),BS->IVtitle[ivVar[0]].Data());
         grName = Form("kindep_%s",BS->IVname[ivVar[0]].Data());
         kindepGr = new TGraphErrors();
         kindepGr->SetName(grName);
         kindepGr->SetTitle(grTitle);
       };

       kindepMap.insert(std::pair<Int_t,TGraphErrors*>(binNum,kindepGr));

     };
   }
   else if(dimensions == 2) {
     for(int b1=0; b1<NB[1]; b1++) {
       for(int b0=0; b0<NB[0]; b0++) {

         A = new Asymmetry(BS, whichModulation, 2, 
           ivVar[0], b0,
           ivVar[1], b1
         );
         asymVec.push_back(A);
         binNum = GetBinNum(b0,b1);
         binMap.insert(std::pair<int,int>(binNum,bcnt));
         asymMap.insert(std::pair<Int_t, Asymmetry*>(binNum,A));
         bcnt++;

         if(b0==0) {
           grTitle = Form("%s asymmetry vs. %s :: %s",
             A->ModulationTitle.Data(),BS->IVtitle[ivVar[0]].Data(),
             (BS->GetBoundStr(ivVar[1],b1)).Data());
           grName = Form("kindep_%s_bin_%s%d",BS->IVname[ivVar[0]].Data(),
             BS->IVname[ivVar[1]].Data(),b1);
           kindepGr = new TGraphErrors();
           kindepGr->SetName(grName);
           kindepGr->SetTitle(grTitle);
         };

         kindepMap.insert(std::pair<Int_t,TGraphErrors*>(binNum,kindepGr));

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
           asymVec.push_back(A);
           binNum = GetBinNum(b0,b1,b2);
           binMap.insert(std::pair<int,int>(binNum,bcnt));
           asymMap.insert(std::pair<Int_t, Asymmetry*>(binNum,A));
           bcnt++;

           if(b0==0) {
             grTitle = Form("%s asymmetry vs. %s :: %s, %s",
               A->ModulationTitle.Data(),BS->IVtitle[ivVar[0]].Data(),
               (BS->GetBoundStr(ivVar[1],b1)).Data(),
               (BS->GetBoundStr(ivVar[2],b2)).Data());
             grName = Form("kindep_%s_bin_%s%d_%s%d",BS->IVname[ivVar[0]].Data(),
               BS->IVname[ivVar[1]].Data(),b1,BS->IVname[ivVar[2]].Data(),b2);
             kindepGr = new TGraphErrors();
             kindepGr->SetName(grName);
             kindepGr->SetTitle(grTitle);
           };

           kindepMap.insert(std::pair<Int_t,TGraphErrors*>(binNum,kindepGr));

         };
       };
     };
   };



   // EVENT LOOP -------------------------------------------
   printf("begin loop through %lld events...\n",ev->ENT);
   Bool_t filled;
   for(int i=0; i<ev->ENT; i++) {

     ev->GetEvent(i);

     if(ev->cutDihadron && ev->cutQ2 && ev->cutW && ev->cutY) {
       
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

         filled = A->FillPlots(); 

         //if(filled && A->debug) ev->PrintEvent();
       };

     };
   };
   // end event loop -------------------------------------------



   // compute asymmetries
   Float_t asymValue,asymError;
   Float_t kinValue,kinError;
   TF1 * fitFunc;
   printf("--- calculate asymmetries\n");
   for(std::vector<Asymmetry*>::iterator it = asymVec.begin(); 
     it!=asymVec.end(); ++it
   ) {

     A = *it;
     A->CalculateAsymmetries();

     fitFunc = A->asymGr->GetFunction("pol1");

     if(fitFunc!=NULL) {
       
       // asymmetry value
       asymValue = fitFunc->GetParameter(1);

       // asymmetry statistical uncertainty
       asymError = fitFunc->GetParError(1);

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

       binNum = GetBinNum(A->B[0], A->B[1], A->B[2]);
       kindepGr = kindepMap.at(binNum);
       kindepGr->SetPoint(A->B[0],kinValue,asymValue);
       kindepGr->SetPointError(A->B[0],kinError,asymError);
     };

   };



   // -- instantiate canvases
   TString canvName = "kindepCanv";
   for(int d=0; d<dimensions; d++) {
     if(d==1) canvName = Form("%s_bins_%s",canvName.Data(),(BS->IVname[ivVar[d]]).Data());
     else canvName = Form("%s_%s",canvName.Data(),(BS->IVname[ivVar[d]]).Data());
   };
   Int_t canvX,canvY,divX,divY;
   Int_t canvSize = 800;
   switch(dimensions) {
     case 1:
       canvX=canvSize; canvY=canvSize;
       divX=1; divY=1;
       break;
     case 2:
       canvX=NB[1]*canvSize; canvY=canvSize;
       divX=NB[1]; divY=1;
       break;
     case 3:
       canvX=3*canvSize; canvY=3*canvSize;
       divX=NB[1]; divY=NB[2];
       break;
   };
   TCanvas * kindepCanv = new TCanvas(canvName,canvName,canvX,canvY);
   kindepCanv->Divide(divX,divY);


   // -- zero line
   //TLine * zeroLine;
   //Float_t drawMin,drawMax;
   /*
   Float_t drawMean,drawRMS;
   switch(dimensions) {
     case 1:
       binNum = GetBinNum(0);
       kindepGr = kindepMap.at(binNum);
       drawMean = (asymMap.at(binNum))->wdist1->
   = new TLine(BS->minIV[ivVar[0]],0,BS->maxIV[ivVar[0]],0);
   zeroLine->SetLineColor(kBlack);
   zeroLine->SetLineWidth(1.5);
   zeroLine->SetLineStyle(kDashed);
   */


   // -- add objects to canvases
   if(dimensions==1) {
     kindepCanv->cd();
     binNum = GetBinNum(0);
     kindepGr = kindepMap.at(binNum);
     DrawKinDepGraph(kindepGr,BS,ivVar[0]);
   }
   else if(dimensions==2) {
     for(int b1=0; b1<NB[1]; b1++) {
       kindepCanv->cd(b1+1);
       binNum = GetBinNum(0,b1);
       kindepGr = kindepMap.at(binNum);
       DrawKinDepGraph(kindepGr,BS,ivVar[0]);
     };
   }
   else if(dimensions==3) {
     for(int b1=0; b1<NB[1]; b1++) {
       for(int b2=0; b2<NB[2]; b2++) {
         kindepCanv->cd(b1*NB[1]+b2+1);
         binNum = GetBinNum(0,b1,b2);
         kindepGr = kindepMap.at(binNum);
         DrawKinDepGraph(kindepGr,BS,ivVar[0]);
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

     A->modDist->Write();
     for(Int_t m=0; m<Asymmetry::nModBins; m++) A->modBinDist[m]->Write();
     if(dimensions==1) A->IVvsModDist->Write();
     for(Int_t s=0; s<nSpin; s++) A->aziDist[s]->Write();
   };

   // -- asymmetries and kindep graphs
   printf("--- write kinematic-dependent asymmetries\n");
   for(std::vector<Asymmetry*>::iterator it = asymVec.begin(); 
     it!=asymVec.end(); ++it
   ) {
     A = *it;
     A->PrintSettings();

     // first write out the asymmetry vs. modulation graphs
     A->asymGr->Write();

     // then write out the kindepGr *after* writing out all the
     // relevant asymmetry vs. modulation graphs
     if(A->B[0] + 1 == NB[0]) {
       binNum = GetBinNum(A->B[0], A->B[1], A->B[2]);
       kindepGr = kindepMap.at(binNum);
       kindepGr->Write();
     };
   };

   kindepCanv->Write();
     

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
