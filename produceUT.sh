#!/bin/bash

buildSpinroot.exe -f outroot/pythia_5x41_smear.root -i1
catSpinroot.exe


asymFit.exe 921 0 | tee fit_test1.log
mv spinroot/asym_921.root spinroot/asym_test1.root

asymFit.exe 801 0 | tee fit_test2.log
mv spinroot/asym_801.root spinroot/asym_test2.root

asymFit.exe 803 0 | tee fit_test3.log
mv spinroot/asym_803.root spinroot/asym_test3.root

mv fit_test*.log spinroot/
