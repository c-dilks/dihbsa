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

  canvSize = 400;


  // 1D canvases
  canvName = "canv1";
  canv1 = new TCanvas(canvName,canvName,
    N*canvSize, canvSize
  );
  canv1->Divide(1,N);
  

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

      canvName = Form("canv_%s_bin_%s",
        A->IVname[v0].Data(),
        A->IVname[v1].Data()
      );
      canv2[v0][v1] = new TCanvas(canvName,canvName,
        NB[v1]*canvSize, canvSize
      );

      for(int b=0; b<NB[v1]; b++) {
        grName = Form("asymGr_%s_bin_%s%d",
          A->IVname[v0],A->IVname[v1],
          b
        );
        grTitle = Form("asymmetry vs. %s :: %s#in[%.2f,%.2f)",
          A->IVname[v0],A->IVname[v1],
          A->bound[v1][b],A->bound[v1][b+1]
        );
        asymGr2[v0][v1][b] = new TGraphErrors();
        asymGr2[v0][v1][b]->SetName(grName);
        asymGr2[v0][v1][b]->SetTitle(grTitle);
      };

    };
  };
  
  // 3D canvases
  for(int v1=0; v1<N; v1++) {
    for(int v2=0; v2<N; v2++) {
      for(int v3=0; v3<N; v3++) {

        canvName = Form("canv_%s_bin_%s_%s",
          A->IVname[v0].Data(),
          A->IVname[v1].Data(),
          A->IVname[v2].Data()
        );
        canv3[v0][v1][v2] = new TCanvas(canvName,canvName,
          3*canvSize, 3*canvSize
        );

      };
    };
  };

  FillCanvases();
        
};


void KinDep::FillCanvases() {
  
  // fill 1D canvases
  for(int v=0; v<N; v++) {
    canv1->cd(v+1);

  };
};




KinDep::~KinDep() {
};

