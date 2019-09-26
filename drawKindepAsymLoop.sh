#!/bin/bash

# loop through some spinFinal*.root files and print
# asymmetry canvases to png files

pushd forCompareKFweights

cat > tmp.C << 'EOF'
void tmp(TString fn) {
  TFile * f = new TFile(fn,"READ");
  //TCanvas * c = (TCanvas*) f->Get("RF_A0_kindepCanv_sinPhiR_M");
  TCanvas * c = (TCanvas*) f->Get("RF_A0_kindepCanv_sinPhiHRweight_M");
  c->Print(TString(fn+".png"),"png");
}
EOF

#for file in e*.root; do
for file in g*.root; do
  root -b -q tmp.C'("'$file'")'
done

rm tmp.C

popd
