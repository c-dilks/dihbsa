R__LOAD_LIBRARY(DihBsa)

#include "Constants.h"
#include "EventTree.h"
#include "Tools.h"

void LookForPi0s(TString dir="outroot.dnp2018.some") {

  TString files = dir + "/*.root";
  TFile * outfile = new TFile("pi0.root","RECREATE");
  Int_t whichPair = EncodePairType(kPip,kPi0); // must have a diphoton
  EventTree * ev = new EventTree(files,whichPair);

  Bool_t cut;
  Bool_t phiCut;

  const Int_t NBINS = 200;
  const Float_t mMax = 0.7;
  const Float_t ptMax = 1.5;
  const Float_t alphaMax = 0.35;
  TH1D * hMass = new TH1D("hMass","diphoton mass;M",NBINS,0,mMax);
  TH1D * hAlpha = new TH1D("hAlpha","diphoton opening angle;#alpha",NBINS,0,alphaMax);
  TH2D * hMAlpha = new TH2D("hMAlpha",
    "diphoton mass vs. opening angle;#alpha;M",
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

    // require dihadron with kinematics and DIS cuts; 
    // one of the hadrons must be of type pi0
    /*
    - this gives us a diphoton, without any pi0 cuts, paired with any of the
      observable hadrons such that the observable-diphoton system satisfies dihadron
      kinematic and DIS cuts
    - in the case where we have {diphoton, pi+, pi-}, there are 3 ways to make
      dihadrons, 2 of which contain the diphoton; this diphoton will be filled twice
      in the distributions
    */
    if( /*ev->cutDihadronKinematics &&*/ ev->cutDIS && 
        (ev->hadIdx[qA]==kPi0 || ev->hadIdx[qB]==kPi0) )
    {

      phiCut = Tools::PhiFiducialCut(ev->photPhi[0]) && 
               Tools::PhiFiducialCut(ev->photPhi[1]);


      // 0
      cut = true;

      // 1
      //cut = ev->photE[0]>0.5 && ev->photE[1]>0.5;

      // 2
      /*
      cut = ev->photE[0]>0.5 && ev->photE[1]>0.5 &&
            ev->diphAlpha > 0.05 && ev->diphAlpha < 0.2;
            */

      // 3
      /*
      cut = ev->photE[0]>0.5 && ev->photE[1]>0.5 &&
            ev->diphAlpha > 0.05 && ev->diphAlpha < 0.2 &&
            ev->diphPt > 0.15;
            */

      // 4
      /*
      cut = ev->photE[0]>0.5 && ev->photE[1]>0.5 &&
            ev->diphAlpha > 0.05 && ev->diphAlpha < 0.2 &&
            ev->diphPt > 0.15 &&
            ev->diphZ > 0.1 && ev->diphZ < 0.6;
            */
      // 5 -- tighten alpha cut and add phi fiducial cut on each photon
      /*
      cut = ev->photE[0]>0.5 && ev->photE[1]>0.5 &&
            ev->diphAlpha > 0.05 && ev->diphAlpha < 0.2 &&
            ev->diphPt > 0.15 &&
            ev->diphZ > 0.1 && ev->diphZ < 0.6 &&
            phiCut;
            */

      // OLD DNP2018 ATTEMPT
      // spans pi0 peak alpha range
      //cut = ev->diphAlpha < 0.3 && ev->diphAlpha > 0.05;
      // now cut out low pT radiative junk
      //cut = ev->diphAlpha < 0.3 && ev->diphAlpha > 0.05 && ev->diphPt>0.15;
      // add on E_hadron > 1 GeV cut (also mass peak is rather wide below 1 GeV)
      //cut = ev->diphAlpha < 0.3 && ev->diphAlpha > 0.05 && ev->diphPt>0.15 &&
            //ev->diphE > 1;
      // standard high-Z cut
      //cut = ev->diphAlpha < 0.3 && ev->diphAlpha > 0.05 && ev->diphPt>0.15 &&
            //ev->diphE > 1 && ev->diphZ < 0.6;
      // make sure both photons are in azimuthal sectors
      //cut = ev->diphAlpha < 0.3 && ev->diphAlpha > 0.05 && ev->diphPt>0.15 &&
            //ev->diphE > 1 && ev->diphZ < 0.6 && phiCut;
      // tighten alpha and energy cuts
      //cut = ev->diphAlpha < 0.25 && ev->diphAlpha > 0.07 && ev->diphPt>0.15 &&
            //ev->diphE > 1.3 && ev->diphZ < 0.6 && phiCut;



      // 5038 CUTS
      //cut = ev->diphE > 2 && ev->diphZ < 0.6; // old cuts
      
      // 0
      //cut = true;
      // 1
      //cut = ev->photE[0]>0.5 && ev->photE[1]>0.5;
      // 2
      /*
      cut = ev->photE[0]>0.5 && ev->photE[1]>0.5 &&
            ev->diphAlpha > 0.04 && ev->diphAlpha < 0.15;
            */

      if(cut) {
        hMass->Fill(ev->diphM);
        hAlpha->Fill(ev->diphAlpha);
        hMAlpha->Fill(ev->diphAlpha,ev->diphM);
        hME->Fill(ev->diphE,ev->diphM);
        hMZ->Fill(ev->diphZ,ev->diphM);
        hMPt->Fill(ev->diphPt,ev->diphM);
        hMEta->Fill(ev->diphEta,ev->diphM);
        hMPhi->Fill(ev->diphPhi,ev->diphM);
        
        photEcorr->Fill(ev->photE[1],ev->photE[0]);
        photPtcorr->Fill(ev->photPt[1],ev->photPt[0]);
        photEtacorr->Fill(ev->photEta[1],ev->photEta[0]);
        photPhicorr->Fill(ev->photPhi[1],ev->photPhi[0]);

        hMphotE->Fill(ev->photE[0],ev->diphM);
        hMphotPt->Fill(ev->photPt[0],ev->diphM);
        hMphotEta->Fill(ev->photEta[0],ev->diphM);
        hMphotPhi->Fill(ev->photPhi[0],ev->diphM);
        photPyPx->Fill(TMath::Cos(ev->photPhi[0])*ev->photPt[0],
                       TMath::Sin(ev->photPhi[0])*ev->photPt[0]);
      };


    };
  };


  hMass->Fit("gaus","","",0.13,0.16);

  hMass->Write();
  hAlpha->Write();
  hMAlpha->Write();
  hME->Write();
  hMZ->Write();
  hMPt->Write();
  hMEta->Write();
  hMPhi->Write();

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

  

