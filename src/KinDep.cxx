#include "KinDep.h"

ClassImp(KinDep)

using namespace std;


KinDep::KinDep(Asymmetry * asym_) {
  printf("Instantiating KinDep...\n");


  A = asym_;
  if(Asymmetry::nIV != 3) {
    fprintf(stderr,"ERROR: KinDep does not know how to process nIV!=3\n");
    return;
  };
  N = Asymmetry::nIV;
  for(int b=0; b<N; b++) NB[b] = A->nBins[b];


  // formatting
  canvSize = 400;

  fitName = "pol1"; // fit function used in Asymmetry fits


  // 1D canvases
  //
  canvName = "canv1";
  canv1 = new TCanvas(canvName,canvName,
    N*canvSize, canvSize
  );
  canv1->Divide(N,1);

  //
  for(int v=0; v<N; v++) {
    grName = Form("asymGr_%s",A->IVname[v].Data());
    grTitle = Form("asymmetry vs. %s",A->IVtitle[v].Data());
    asymGr1[v] = new TGraphErrors();
    asymGr1[v]->SetName(grName);
    asymGr1[v]->SetTitle(grTitle);
  };

  
  // 2D canvases
  for(int v0=0; v0<N; v0++) {
    for(int v1=0; v1<N; v1++) {

      //
      canvName = Form("canv_%s_bin_%s",
        A->IVname[v0].Data(),
        A->IVname[v1].Data()
      );
      canv2[v0][v1] = new TCanvas(canvName,canvName,
        NB[v1]*canvSize, canvSize
      );
      canv2[v0][v1]->Divide(NB[v1],1);

      //
      for(int b=0; b<NB[v1]; b++) {
        grName = Form("asymGr_%s_bin_%s%d",
          A->IVname[v0].Data(),A->IVname[v1].Data(),
          b
        );
        grTitle = Form("asymmetry vs. %s :: %s#in[%.2f,%.2f)",
          A->IVname[v0].Data(),
          A->IVname[v1].Data(),A->bound[v1][b],A->bound[v1][b+1]
        );
        asymGr2[v0][v1][b] = new TGraphErrors();
        asymGr2[v0][v1][b]->SetName(grName);
        asymGr2[v0][v1][b]->SetTitle(grTitle);
      };

    };
  };
  
  // 3D canvases
  for(int v0=0; v0<N; v0++) {
    for(int v1=0; v1<N; v1++) {
      for(int v2=0; v2<N; v2++) {

        //
        canvName = Form("canv_%s_bin_%s_%s",
          A->IVname[v0].Data(),
          A->IVname[v1].Data(),
          A->IVname[v2].Data()
        );
        canv3[v0][v1][v2] = new TCanvas(canvName,canvName,
          3*canvSize, 3*canvSize
        );
        canv3[v0][v1][v2]->Divide(NB[v1],NB[v2]);

        //
        for(int b1=0; b1<NB[v1]; b1++) {
          for(int b2=0; b2<NB[v1]; b2++) {
            grName = Form("asymGr_%s_bin_%s%d_%s%d",
              A->IVname[v0].Data(),
              A->IVname[v1].Data(),b1,
              A->IVname[v2].Data(),b2
            );
            grTitle = Form("asymmetry vs. %s :: %s#in[%.2f,%.2f),  %s#in[%.2f,%.2f)",
              A->IVname[v0].Data(),
              A->IVname[v1].Data(),A->bound[v1][b1],A->bound[v1][b1+1],
              A->IVname[v2].Data(),A->bound[v2][b2],A->bound[v2][b2+1]
            );
            asymGr3[v0][v1][v2][b1][b2] = new TGraphErrors();
            asymGr3[v0][v1][v2][b1][b2]->SetName(grName);
            asymGr3[v0][v1][v2][b1][b2]->SetTitle(grTitle);
          };
        };
      };
    };
  };

  FillAsymGraphs();
  FillCanvases();
        
};


void KinDep::FormatAsymGr(TGraphErrors * g, Int_t ivNum) {

  switch(ivNum) {
    case Asymmetry::vM: g->SetLineColor(kRed); break;
    case Asymmetry::vX: g->SetLineColor(kGreen+1); break;
    case Asymmetry::vZ: g->SetLineColor(kBlue); break;
    default: return;
  };

  g->SetLineWidth(2);

  g->SetMarkerStyle(kFullCircle);
  g->SetMarkerColor(kBlack);
  g->SetMarkerSize(1.3);

  g->GetYaxis()->SetRangeUser(-0.05,0.1);
};

  

void KinDep::FillAsymGraphs() {

  // fill 1D graphs
  for(int v=0; v<N; v++) {
    for(int p=0; p<NB[v]; p++) {

      fitFunc = A->asym1[v][p]->GetFunction(fitName);

      asymValue = fitFunc->GetParameter(1);
      asymError = fitFunc->GetParError(1);

      kinValue = A->wDist1[v][p]->GetMean();
      kinError = A->wDist1[v][p]->GetRMS();

      asymGr1[v]->SetPoint(p,kinValue,asymValue);
      asymGr1[v]->SetPointError(p,kinError,asymError);
    };
  };


  // fill 2D graphs
  for(int v0=0; v0<N; v0++) {
    for(int v1=0; v1<N; v1++) {
      for(int b1=0; b1<NB[v1]; b1++) {

        for(int p=0; p<NB[v0]; p++) {

          fitFunc = A->asym2[v0][v1][p][b1]->GetFunction(fitName);

          asymValue = fitFunc->GetParameter(1);
          asymError = fitFunc->GetParError(1);

          kinValue = A->wDist2[v0][v1][p][b1]->GetMean(1);
          kinError = A->wDist2[v0][v1][p][b1]->GetRMS(1);

          asymGr2[v0][v1][b1]->SetPoint(p,kinValue,asymValue);
          asymGr2[v0][v1][b1]->SetPointError(p,kinError,asymError);
        };
      };
    };
  };


  // fill 3D graphs
  /*
  for(int v0=0; v0<N; v0++) {
    for(int v1=0; v1<N; v1++) {
      for(int v2=0; v2<N; v2++) {
        for(int b1=0; b1<NB[v1]; b1++) {
          for(int b2=0; b2<NB[v2]; b2++) {

            for(int p=0; p<NB[v0]; p++) {

              printf("aqui %d %d %d %d %d %d\n",v0,v1,v2,b1,b2,p);

              // TODO -- check this assigment
              switch(v0) {
                case Asymmetry::vM: 
                  asymGrCurr = A->asym3[p][b1][b2];
                  wDistCurr = A->wDist3[p][b1][b2];
                  break;
                case Asymmetry::vX: 
                  asymGrCurr = A->asym3[b1][p][b2];
                  wDistCurr = A->wDist3[b1][p][b2];
                  break;
                case Asymmetry::vZ: 
                  asymGrCurr = A->asym3[b1][b2][p];
                  wDistCurr = A->wDist3[b1][b2][p];
                  break;
              };



              printf("aqui switch %p\n",(void*)asymGrCurr);
              fitFunc = asymGrCurr->GetFunction(fitName);
              printf("aqui switch %p\n",(void*)fitFunc);

              asymValue = fitFunc->GetParameter(1);
              asymError = fitFunc->GetParError(1);

              kinValue = wDistCurr->GetMean(v0+1);
              kinError = wDistCurr->GetRMS(v0+1);

              printf("aqui vals\n");

              asymGr3[v0][v1][v2][b1][b2]->SetPoint(p,kinValue,asymValue);
              asymGr3[v0][v1][v2][b1][b2]->SetPointError(p,kinError,asymError);
              printf("aqui fills\n");
            };
          };
        };
      };
    };
  };
  */

};




void KinDep::FillCanvases() {
  
  // fill 1D canvases
  for(int v=0; v<N; v++) {
    canv1->cd(v+1);
    asymGr1[v]->Draw("APE");
    FormatAsymGr(asymGr1[v],v);
    asymGr1[v]->Draw("APE");
  };


  // fill 2D canvases
  for(int v0=0; v0<N; v0++) {
    for(int v1=0; v1<N; v1++) {
      for(int b1=0; b1<NB[v1]; b1++) {
        canv2[v0][v1]->cd(b1+1);
        asymGr2[v0][v1][b1]->Draw("APE");
        FormatAsymGr(asymGr2[v0][v1][b1],v0);
        asymGr2[v0][v1][b1]->Draw("APE");
      };
    };
  };

  // fill 3D canvases
  /*
  for(int v0=0; v0<N; v0++) {
    for(int v1=0; v1<N; v1++) {
      for(int v2=0; v2<N; v2++) {
        for(int b1=0; b1<NB[v1]; b1++) {
          for(int b2=0; b2<NB[v2]; b2++) {

            // fill 2d grid of canvases, such that
            // -- bottom-left corner is for ( b1=0, b2=0 )
            // -- top-left corner is for ( b1=0, b2=NB[v2]-1 )
            canv3[v0][v1][v2]->cd(
              ( NB[v2] - 1 - b2 ) * NB[v1] + b1 + 1
            );
            asymGr3[v0][v1][v2][b1][b2]->Draw("APE");
            FormatAsymGr(asymGr3[v0][v1][v2][b1][b2],v0);
            asymGr3[v0][v1][v2][b1][b2]->Draw("APE");
          };
        };
      };
    };
  };
  */


};


// to be called AFTER Asymmetry::Write()
void KinDep::Write(TFile * f) {
  f->cd();

  // 1D write
  f->cd("/1D");
  canv1->Write();

  // 2D write
  f->cd("/"); f->mkdir("2D/canv"); f->cd("/2D/canv");
  for(int v0=0; v0<N; v0++) {
    for(int v1=0; v1<N; v1++) {
      if(v0!=v1) canv2[v0][v1]->Write();
    };
  };

  // 3D write
  /*
  f->cd("/"); f->mkdir("3D/canv"); f->cd("/3D/canv");
  for(int v0=0; v0<N; v0++) {
    for(int v1=0; v1<N; v1++) {
      for(int v2=0; v2<N; v2++) {
        canv3[v0][v1][v2]->Write();
      };
    };
  };
  */
};



KinDep::~KinDep() {
};

