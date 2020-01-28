// vary Dpararm for evaluating systematic uncertainty on asymmetry from unpolarized DiFF
// pp-wave (unmeasured)

// OPTIONS
// --------------
Int_t fitMode = 51; // Asymmetry.cxx fitMode
Float_t VMfrac = 0.12; // fraction of dihadrons from Vector Mesons (evaluated in MC)
Float_t nsteps = 5; // number of steps of Dparam to vary
// --------------

TString cmdFit,cmdMv;

void exeFit(Float_t d) {
  cmdFit = Form("asymFit.exe %d spinroot %f",fitMode,d);
  cmdMv = Form("mv -v spinroot/asym_%d.root spinroot/Dasym_%d_%f.root",fitMode,fitMode,d);
  system(cmdFit);
  system(cmdMv);
  system("sleep 3");
};

void VaryDparam() {
  // - this is used as  in the positivity bounds

  // positivity bounds: uses D_{1,OO}^{p}/D_{1,OO} set to 4/3 * VMfrac
  VMfrac *= 4.0/3.0;
  Float_t lb = -3.0/2.0 * VMfrac;
  Float_t ub = 3.0 * VMfrac;
  printf("VMfrac = %f\n",VMfrac);
  printf("D_{1,LL} / D_{1,OO} will be varied between %f and %f\n",lb,ub);

  // determine sequence of Dparam values
  Float_t stepsize = (ub-lb)/nsteps;

  // loop over the sequence, executing the fits
  Float_t dval = lb;
  while(dval <= ub) {
    printf("EXECUTE FIT: dval = %f\n",dval);
    system("sleep 1");
    exeFit(dval);
    dval += stepsize;
  };
  printf("EXECUTE DEFAULT FIT\n");
  system("sleep 1");
  exeFit(0); // execute default fit (without Dparam)
};

