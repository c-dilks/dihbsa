#!/bin/bash
# wrapper for streamlining RatioPlot.C

outdir="PperpStudy3/z_between_0.2_and_0.3"
mkdir -p $outdir; rm -r $outdir; mkdir -p $outdir

root -b -q RatioPlot.C'("roots_z_below_0.3/eicPlots.pythia_5x41_smear.ymin_0.03.root","roots_z_below_0.3/eicPlots.pythia_5x41_smear.ymin_0.00.root","hadronA","Pperp_PiPlus_ymin0.03","PperpDistLin")'
root -b -q RatioPlot.C'("roots_z_below_0.3/eicPlots.pythia_5x41_smear.ymin_0.03.root","roots_z_below_0.3/eicPlots.pythia_5x41_smear.ymin_0.00.root","hadronA","Qt_PiPlus_ymin0.03","QtDistLin")'
root -b -q RatioPlot.C'("roots_z_below_0.3/eicPlots.pythia_5x41_smear.ymin_0.03.root","roots_z_below_0.3/eicPlots.pythia_5x41_smear.ymin_0.00.root","hadronA","QtOverQ_PiPlus_ymin0.03","QtOverQdist")'

root -b -q RatioPlot.C'("roots_z_below_0.3/eicPlots.pythia_5x41_smear.ymin_0.05.root","roots_z_below_0.3/eicPlots.pythia_5x41_smear.ymin_0.00.root","hadronA","Pperp_PiPlus_ymin0.05","PperpDistLin")'
root -b -q RatioPlot.C'("roots_z_below_0.3/eicPlots.pythia_5x41_smear.ymin_0.05.root","roots_z_below_0.3/eicPlots.pythia_5x41_smear.ymin_0.00.root","hadronA","Qt_PiPlus_ymin0.05","QtDistLin")'
root -b -q RatioPlot.C'("roots_z_below_0.3/eicPlots.pythia_5x41_smear.ymin_0.05.root","roots_z_below_0.3/eicPlots.pythia_5x41_smear.ymin_0.00.root","hadronA","QtOverQ_PiPlus_ymin0.05","QtOverQdist")'

#mv *PiPlus_ymin*bin{1,4,7,10,13}.pdf $outdir/
mv *PiPlus_ymin*bin*.pdf $outdir/



outdir="PperpStudy3/z_between_0.3_and_1.0"
mkdir -p $outdir; rm -r $outdir; mkdir -p $outdir

root -b -q RatioPlot.C'("roots_z_above_0.3/eicPlots.pythia_5x41_smear.ymin_0.03.root","roots_z_above_0.3/eicPlots.pythia_5x41_smear.ymin_0.00.root","hadronA","Pperp_PiPlus_ymin0.03","PperpDistLin")'
root -b -q RatioPlot.C'("roots_z_above_0.3/eicPlots.pythia_5x41_smear.ymin_0.03.root","roots_z_above_0.3/eicPlots.pythia_5x41_smear.ymin_0.00.root","hadronA","Qt_PiPlus_ymin0.03","QtDistLin")'
root -b -q RatioPlot.C'("roots_z_above_0.3/eicPlots.pythia_5x41_smear.ymin_0.03.root","roots_z_above_0.3/eicPlots.pythia_5x41_smear.ymin_0.00.root","hadronA","QtOverQ_PiPlus_ymin0.03","QtOverQdist")'

root -b -q RatioPlot.C'("roots_z_above_0.3/eicPlots.pythia_5x41_smear.ymin_0.05.root","roots_z_above_0.3/eicPlots.pythia_5x41_smear.ymin_0.00.root","hadronA","Pperp_PiPlus_ymin0.05","PperpDistLin")'
root -b -q RatioPlot.C'("roots_z_above_0.3/eicPlots.pythia_5x41_smear.ymin_0.05.root","roots_z_above_0.3/eicPlots.pythia_5x41_smear.ymin_0.00.root","hadronA","Qt_PiPlus_ymin0.05","QtDistLin")'
root -b -q RatioPlot.C'("roots_z_above_0.3/eicPlots.pythia_5x41_smear.ymin_0.05.root","roots_z_above_0.3/eicPlots.pythia_5x41_smear.ymin_0.00.root","hadronA","QtOverQ_PiPlus_ymin0.05","QtOverQdist")'

#mv *PiPlus_ymin*bin{1,4,7,10,13}.pdf $outdir/
mv *PiPlus_ymin*bin*.pdf $outdir/
