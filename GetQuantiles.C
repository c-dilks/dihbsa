// computes quantiles for x,Mh,Zpair
// nQuant quantiles will be determined

void GetQuantiles(TString plotsFile="plots.root", Int_t nQuant = 5) {
  const Int_t N = nQuant;
  Double_t q[N];
  Double_t p[N];
  int i;
  for(i=0; i<N; i++) p[i] = Double_t(i+1)/N;

  TFile * f = new TFile(plotsFile,"READ");
  TH1D * XDist = (TH1D*) f->Get("XDist");
  TH1D * MhDist = (TH1D*) f->Get("MhDist");
  TH1D * ZpairDist = (TH1D*) f->Get("ZpairDist");

  XDist->GetQuantiles(N,q,p);
  printf("XDist:\n");
  for(i=0; i<N; i++) printf("%f\n",q[i]);
  printf("\n");

  MhDist->GetQuantiles(N,q,p);
  printf("MhDist:\n");
  for(i=0; i<N; i++) printf("%f\n",q[i]);
  printf("\n");

  ZpairDist->GetQuantiles(N,q,p);
  printf("ZpairDist:\n");
  for(i=0; i<N; i++) printf("%f\n",q[i]);
  printf("\n");
};
