R__LOAD_LIBRARY(DihBsa)

#include "Constants.h"

void DiagMatrix() {

  gStyle->SetOptStat(0);

  Int_t dim = 1;
  //TString plotName = "ZpairDist";
  //TString plotName = "XDist";
  //TString plotName = "MhDist";
  //TString plotName = "thetaDist";
  TString plotName = "PhiHRDist";
  //Int_t dim = 2;
  //TString plotName = "PhiHvsPhiR";


  TString canvName = "matrix_" + plotName;
  TCanvas * matrix = new TCanvas(canvName,canvName,1000,1000);
  matrix->Divide(nObservables,nObservables);

  TFile * infile[nObservables][nObservables];
  TH1D * dist1[nObservables][nObservables];
  TH2D * dist2[nObservables][nObservables];


  TString pairN,pairT;
  TString fileN;
  TString titleTmp;
  Int_t pa,pb;
  Int_t pad;
  for(int o1=0; o1<nObservables; o1++) {
    for(int o2=o1; o2<nObservables; o2++) {
      IterPair(o1,o2,pa,pb);
      pairN = PairName(pa,pb);
      pairT = PairTitle(pa,pb);

      pad = (nObservables-1-o1)*nObservables + o2 + 1;

      printf("obs(%d,%d) = (%s,%s) = \"%s\"   pad=%d\n",
        o1,o2,
        ObsName(o1).Data(),ObsName(o2).Data(),pairN.Data(),
        pad);

      fileN = "diagset/plots." + pairN + ".root";
      infile[o1][o2] = new TFile(fileN,"READ");

      matrix->cd(pad);
      switch(dim) {
        case 1:
          matrix->GetPad(pad)->SetGrid(1,0);
          dist1[o1][o2] = (TH1D*) infile[o1][o2]->Get(plotName);
          titleTmp = dist1[o1][o2]->GetTitle();
          titleTmp = pairT + " " + titleTmp;
          dist1[o1][o2]->SetTitle(titleTmp);
          dist1[o1][o2]->SetLineWidth(3);
          dist1[o1][o2]->GetXaxis()->SetRangeUser(-PI,PI);
          dist1[o1][o2]->Draw();
          break;
        case 2:
          matrix->GetPad(pad)->SetGrid(1,1);
          dist2[o1][o2] = (TH2D*) infile[o1][o2]->Get(plotName);
          titleTmp = dist2[o1][o2]->GetTitle();
          titleTmp = pairT + " " + titleTmp;
          dist2[o1][o2]->SetTitle(titleTmp);
          //dist2[o1][o2]->GetXaxis()->SetRangeUser(-PI,PI);
          //dist2[o1][o2]->GetYaxis()->SetRangeUser(-PI,PI);
          dist2[o1][o2]->Draw("colz");
          break;
      };
    };
  };
};

