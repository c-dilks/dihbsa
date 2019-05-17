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
#include "TCanvas.h"

// DihBsa
#include "Constants.h"
#include "DIS.h"
#include "Trajectory.h"
#include "Dihadron.h"
#include "EventTree.h"
#include "Binning.h"
#include "Asymmetry.h"

Int_t GetBinNum(Int_t bin0, Int_t bin1=-1, Int_t bin2=-1);


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
   /*
   std::map<Int_t, Asymmetry*> asymMap; // map of asymVec entry number (iterator) 
                                        // and Asymmetry instance (this may be slower
                                        // than asymVec, but protects against iterator
                                        // invalidation)
                                        */
   std::map<int,int> binMap; // 3-digit bin number -> asym iterator
                                 // -- use GetBinNum(...) to convert a set of bin
                                 //    numbers to the 3-digit bin number; then use this
                                 //    map to get the index of the corresponding
                                 //    Asymmetry pointer in asymVec

   Asymmetry * A;

   Int_t binNum;
   Int_t bcnt = 0;
   if(dimensions == 1) {
     for(int b=0; b<NB[0]; b++) {
       A = new Asymmetry(BS, whichModulation, 1, ivVar[0], b);
       asymVec.push_back(A);
       binNum = GetBinNum(b);
       binMap.insert(std::pair<int,int>(binNum,bcnt));
       //asymMap.insert(std::pair(binNum,A));
       bcnt++;

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
         //asymMap.insert(std::pair(binNum,A));
         bcnt++;
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
           //asymMap.insert(std::pair(binNum,A));
           bcnt++;
         };
       };
     };
   };



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



   for(std::vector<Asymmetry*>::iterator it = asymVec.begin(); 
     it!=asymVec.end(); ++it
   ) {

     A = *it;
     A->CalculateAsymmetries();

     switch(dimensions) {
       case 1: A->ivDist1->Write(); break;
       case 2: A->ivDist2->Write(); break;
       case 3: A->ivDist3->Write(); break;
     };

     A->modDist->Write();
     for(Int_t m=0; m<Asymmetry::nModBins; m++) A->modBinDist[m]->Write();
     if(dimensions==1) A->IVvsModDist->Write();
     for(Int_t s=0; s<nSpin; s++) A->aziDist[s]->Write();
     A->asymGr->Write();
   };

   outfile->Close();

   printf("end %s\n",argv[0]);

};


Int_t GetBinNum(Int_t bin0, Int_t bin1, Int_t bin2) {
  Int_t retval = bin0;
  if(bin1>=0) retval += 10 * bin1;
  if(bin2>=0) retval += 100 * bin2;
  return retval;
};


