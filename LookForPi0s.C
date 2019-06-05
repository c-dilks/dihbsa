R__LOAD_LIBRARY(DihBsa)

#include "Constants.h"
#include "EventTree.h"
#include "Tools.h"

void LookForPi0s(TString dir="outroot.dnp2018.some", Int_t whichPair=pairP0) {

  TString files = dir + "/*.root";
  TFile * outfile = new TFile("pi0.root","RECREATE");
  EventTree * ev = new EventTree(files,pairP0);

  Bool_t cut;
  Bool_t phiCut;

  const Int_t NBINS = 200;
  const Float_t mMax = 0.7;
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
    NBINS,0,0.8,NBINS,0,mMax);
  TH2D * hMEta = new TH2D("hMEta",
    "diphoton mass vs. pseudorapidity;#eta;M",
    NBINS,0,7,NBINS,0,mMax);
  TH2D * hMPhi = new TH2D("hMPhi",
    "diphoton mass vs. azimuth;#phi;M",
    NBINS,-PIe,PIe,NBINS,0,mMax);


  for(int i=0; i<ev->ENT; i++) {
    ev->GetEvent(i);


    phiCut = Tools::PhiFiducialCut(ev->diphPhi);


    //cut = true;
    
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


    // make sure pion is in one of the sectors (actually this cut ought to be done
    // on each of the photons...)
    //cut = ev->diphAlpha < 0.3 && ev->diphAlpha > 0.05 && ev->diphPt>0.15 &&
          //ev->diphE > 1 && ev->diphZ < 0.6 && phiCut;

    // tighten alpha and energy cuts
    cut = ev->diphAlpha < 0.25 && ev->diphAlpha > 0.07 && ev->diphPt>0.15 &&
          ev->diphE > 1.3 && ev->diphZ < 0.6 && phiCut;




          // CUT ON PHOTON ENERGIES !!!!!!!!!!!!!!!!!!!!!!!
          // FOCAL photons
          // THETA (com frame dihadron angle)
          // Kaons
          // MC
          // bg corrections
          //

    

    //cut = ev->diphE > 2 && ev->diphZ < 0.6; // 5038 skim file cuts

    if(cut) {
      hMass->Fill(ev->diphM);
      hAlpha->Fill(ev->diphAlpha);
      hMAlpha->Fill(ev->diphAlpha,ev->diphM);
      hME->Fill(ev->diphE,ev->diphM);
      hMZ->Fill(ev->diphZ,ev->diphM);
      hMPt->Fill(ev->diphPt,ev->diphM);
      hMEta->Fill(ev->diphEta,ev->diphM);
      hMPhi->Fill(ev->diphPhi,ev->diphM);
    };
  };


  hMass->Fit("gaus","","",0.105,0.16);

  hMass->Write();
  hAlpha->Write();
  hMAlpha->Write();
  hME->Write();
  hMZ->Write();
  hMPt->Write();
  hMEta->Write();
  hMPhi->Write();

};

  

