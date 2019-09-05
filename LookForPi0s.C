R__LOAD_LIBRARY(DihBsa)

#include "Constants.h"
#include "EventTree.h"
#include "Tools.h"
#include "TLorentzVector.h"
#include "TString.h"
#include "TMath.h"

void LookForPi0s(TString dir="outroot.fall18.some") {

  TString files = dir + "/*.root";
  TFile * outfile = new TFile("pi0.root","RECREATE");

  // instantiate EventTree; it doesn't matter what whichPair is, because the cuts that
  // we use from EventTree do not involve selecting pairType (we want to look at all the
  // diphotons "inclusively")
  Int_t whichPair = EncodePairType(kPip,kDiph);
  EventTree * ev = new EventTree(files,whichPair);

  Bool_t cut;
  Bool_t phiCut;
  Float_t angEle[2];
  TLorentzVector photMom[2];

  const Int_t NBINS = 200;
  const Float_t mMax = 0.7;
  const Float_t ptMax = 1.5;
  const Float_t alphaMax = 50;
  TH1D * hMass = new TH1D("hMass","diphoton mass;M",NBINS,0,mMax);
  TH1D * hAlpha = new TH1D("hAlpha","diphoton opening angle [deg];#alpha",
    NBINS,0,alphaMax);
  TH2D * hMAlpha = new TH2D("hMAlpha",
    "diphoton mass vs. opening angle [deg];#alpha;M",
    NBINS,0,alphaMax,NBINS,0,mMax);
  TH2D * hME = new TH2D("hME",
    "diphoton mass vs. energy;E;M",
    NBINS,0,7,NBINS,0,mMax);
  TH2D * hMZ = new TH2D("hMZ",
    "diphoton mass vs. energy sharing (Z=|E_{1}-E_{2}|/E);Z;M",
    NBINS,0,1,NBINS,0,mMax);
  TH2D * hMPt = new TH2D("hMPt",
    "diphoton mass vs. transverse momentum;p_{T};M",
    NBINS,0,ptMax,NBINS,0,mMax);
  TH2D * hMEta = new TH2D("hMEta",
    "diphoton mass vs. pseudorapidity;#eta;M",
    NBINS,0,7,NBINS,0,mMax);
  TH2D * hMPhi = new TH2D("hMPhi",
    "diphoton mass vs. azimuth;#phi;M",
    NBINS,-PIe,PIe,NBINS,0,mMax);

  TH2D * hMAngEle = new TH2D("hMangle",
    "diphoton mass vs. #theta_{e#gamma} [deg];#theta_{e#gamma};M",
    NBINS,0,30,NBINS,0,mMax);

  TH2D * photEcorr = new TH2D("photEcorr","photon E correlation",
    NBINS,0,7,NBINS,0,7);
  TH2D * photPtcorr = new TH2D("photPtcorr","photon p_{T} correlation",
    NBINS,0,ptMax,NBINS,0,ptMax);
  TH2D * photEtacorr = new TH2D("photEtacorr","photon #eta correlation",
    NBINS,0,7,NBINS,0,7);
  TH2D * photPhicorr = new TH2D("photPhicorr","photon #phi correlation",
    NBINS,-PIe,PIe,NBINS,-PIe,PIe);

  TH2D * hMphotE = new TH2D("hMphotE","M_{#gamma#gamma} vs. photon1 E",
    NBINS,0,7,NBINS,0,mMax);
  TH2D * hMphotPt = new TH2D("hMphotPt","M_{#gamma#gamma} vs. photon1 p_{T}",
    NBINS,0,ptMax,NBINS,0,mMax);
  TH2D * hMphotEta = new TH2D("hMphotEta","M_{#gamma#gamma} vs. photon1 #eta",
    NBINS,0,7,NBINS,0,mMax);
  TH2D * hMphotPhi = new TH2D("hMphotPhi","M_{#gamma#gamma} vs. photon1 #phi",
    NBINS,-PIe,PIe,NBINS,0,mMax);

  TH2D * photPyPx = new TH2D("photPyPx","photon1 p_{y} vs. p_{x}",
    NBINS,-ptMax,ptMax,NBINS,-ptMax,ptMax);


  for(int i=0; i<ev->ENT; i++) {
    ev->GetEvent(i);

    // check basic cuts:
    // - cutDIS ensures that x,Q,W are in proper range
    // - cutDihadronKinematics checks Z,Mmiss,xF,hadP (but *not* pairType)
    if( ev->cutDihadronKinematics && ev->cutDIS ) {

      // loop over the diphotons in the hadron
      // (if there are none, diphCnt==0)
      for(int h=0; h<ev->diphCnt; h++) {

        // check if photons are in fiducial phi range
        phiCut = Tools::PhiFiducialCut(ev->diphPhotPhi[h][0]) && 
                 Tools::PhiFiducialCut(ev->diphPhotPhi[h][1]);

        // get angle between photon and electron
        for(int p=0; p<2; p++) {
          photMom[p].SetPtEtaPhiE(
            ev->diphPhotPt[h][p], ev->diphPhotEta[h][p],
            ev->diphPhotPhi[h][p], ev->diphPhotE[h][p] );
          angEle[p] = ev->GetAngleWrtElectron(photMom[p]);
          angEle[p] *= TMath::RadToDeg();
        };


        // 5038 CUTS
        //cut = ev->diphE[h] > 2 && ev->diphZ[h] < 0.6; // old cuts
        
        // 0
        //cut = true;
        // 1 // angEle cuts
        //cut = angEle[0]>5 && angEle[1]>5;
        // 2 // photE cuts
        //cut = ev->diphPhotE[h][0]>0.5 && ev->diphPhotE[h][1]>0.5;
        // 3 // combine angEle and photE cuts
        cut = angEle[0]>5 && angEle[1]>5 &&
              ev->diphPhotE[h][0]>0.5 && ev->diphPhotE[h][1]>0.5;
        // 4 // tighten photE cuts
        //cut = angEle[0]>5 && angEle[1]>5 &&
              //ev->diphPhotE[h][0]>1 && ev->diphPhotE[h][1]>1;
        // 5 // tighten angEle cuts
        //cut = angEle[0]>8 && angEle[1]>8 &&
              //ev->diphPhotE[h][0]>1 && ev->diphPhotE[h][1]>1;




        // TUNING FOR DNP2018 DATASET
        // 0
        //cut = true;

        // 1
        //cut = ev->diphPhotE[h][0]>0.5 && ev->diphPhotE[h][1]>0.5;

        // 2
        /*
        cut = ev->diphPhotE[h][0]>0.5 && ev->diphPhotE[h][1]>0.5 &&
              ev->diphAlpha[h] > 0.05 && ev->diphAlpha[h] < 0.2;
              */

        // 3
        /*
        cut = ev->diphPhotE[h][0]>0.5 && ev->diphPhotE[h][1]>0.5 &&
              ev->diphAlpha[h] > 0.05 && ev->diphAlpha[h] < 0.2 &&
              ev->diphPt[h] > 0.15;
              */

        // 4
        /*
        cut = ev->diphPhotE[h][0]>0.5 && ev->diphPhotE[h][1]>0.5 &&
              ev->diphAlpha[h] > 0.05 && ev->diphAlpha[h] < 0.2 &&
              ev->diphPt[h] > 0.15 &&
              ev->diphZ[h] > 0.1 && ev->diphZ[h] < 0.6;
              */
        // 5 -- tighten alpha cut and add phi fiducial cut on each photon
        /*
        cut = ev->diphPhotE[h][0]>0.5 && ev->diphPhotE[h][1]>0.5 &&
              ev->diphAlpha[h] > 0.05 && ev->diphAlpha[h] < 0.2 &&
              ev->diphPt[h] > 0.15 &&
              ev->diphZ[h] > 0.1 && ev->diphZ[h] < 0.6 &&
              phiCut;
              */




        if(cut) {
          hMass->Fill(ev->diphM[h]);
          hAlpha->Fill(ev->diphAlpha[h]*TMath::RadToDeg());
          hMAlpha->Fill(ev->diphAlpha[h]*TMath::RadToDeg(),ev->diphM[h]);
          hME->Fill(ev->diphE[h],ev->diphM[h]);
          hMZ->Fill(ev->diphZ[h],ev->diphM[h]);
          hMPt->Fill(ev->diphPt[h],ev->diphM[h]);
          hMEta->Fill(ev->diphEta[h],ev->diphM[h]);
          hMPhi->Fill(ev->diphPhi[h],ev->diphM[h]);
          for(int p=0; p<2; p++) hMAngEle->Fill(angEle[p],ev->diphM[h]);

          photEcorr->Fill(ev->diphPhotE[h][1],ev->diphPhotE[h][0]);
          photPtcorr->Fill(ev->diphPhotPt[h][1],ev->diphPhotPt[h][0]);
          photEtacorr->Fill(ev->diphPhotEta[h][1],ev->diphPhotEta[h][0]);
          photPhicorr->Fill(ev->diphPhotPhi[h][1],ev->diphPhotPhi[h][0]);

          hMphotE->Fill(ev->diphPhotE[h][0],ev->diphM[h]);
          hMphotPt->Fill(ev->diphPhotPt[h][0],ev->diphM[h]);
          hMphotEta->Fill(ev->diphPhotEta[h][0],ev->diphM[h]);
          hMphotPhi->Fill(ev->diphPhotPhi[h][0],ev->diphM[h]);
          photPyPx->Fill(TMath::Cos(ev->diphPhotPhi[h][0])*ev->diphPhotPt[h][0],
                         TMath::Sin(ev->diphPhotPhi[h][0])*ev->diphPhotPt[h][0]);
        };


      };
    };
  };


  hMass->Fit("gaus","","",0.110,0.155);

  hMass->Write();
  hAlpha->Write();
  hMAlpha->Write();
  hME->Write();
  hMZ->Write();
  hMPt->Write();
  hMEta->Write();
  hMPhi->Write();
  hMAngEle->Write();

  photEcorr->Write();
  photPtcorr->Write();
  photEtacorr->Write();
  photPhicorr->Write();

  hMphotE->Write();
  hMphotPt->Write();
  hMphotEta->Write();
  hMphotPhi->Write();

  photPyPx->Write();

};

  

