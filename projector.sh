#!/bin/bash
#calls Projector.C for various cases

# projection for hL(x)
root -b -q Projector.C\
'("spinroot.proj.vsX/asym_4.root","xdep",30,1,0.85,0.2)'
mv kindepMA_A2_X.proj.pdf projectionX.pdf
root -b -q Projector.C\
'("spinroot.proj.vsX/asym_4.root","xdep2",30,1,0.85,0.2)'
mv kindepMA_A2_X.proj.pdf projectionX2.pdf

# projection for G1perp
root -b -q Projector.C\
'("spinroot.proj.vsM/asym_4.root","spec",30,1,0.85,0.2)'
mv kindepMA_A1_M.proj.pdf projectionM_30days.pdf
root -b -q Projector.C\
'("spinroot.proj.vsM/asym_4.root","spec",90,1,0.85,0.2)'
mv kindepMA_A1_M.proj.pdf projectionM_90days.pdf

# projection for DSA
root -b -q Projector.C\
'("spinroot.proj.vsZ/asym_1000.root","zero",30,0.85,0.85,0.2)'
mv kindepMA_A1_Z.proj.pdf projectionZ.pdf
#root -b -q Projector.C\
#'("spinroot.proj.vsZ/asym_1001.root","zero",30,0.85,0.85,0.2)'
#mv kindepMA_A0_Z.proj.pdf projectionPW0.pdf
#mv kindepMA_A1_Z.proj.pdf projectionPW1.pdf
#mv kindepMA_A2_Z.proj.pdf projectionPW2.pdf
