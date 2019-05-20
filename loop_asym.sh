#!/bin/bash

dir="outroot.dnp2018"
modulation=0
phiR=3

if [ $# -ge 1 ]; then dir=$1; fi
if [ $# -ge 2 ]; then modulation=$2; fi
if [ $# -ge 3 ]; then phiR=$3; fi

exe="asym.exe"
outdir="spinout"
allCombos=0
argfile="args.dat"
> $argfile

# total number of IVs
nIV=4


mkdir -p $outdir


# 1D
for ((iv=0; iv<$nIV; iv++)); do
  echo $iv
  echo $dir $modulation 1 $iv $phiR 1 >> $argfile
done


# 2D
for ((iv1=0; iv1<$nIV; iv1++)); do
  if [ $allCombos -eq 1 ]; then start=$[$iv1+1];
  else start=0; fi
  for ((iv2=$start; iv2<$nIV; iv2++)); do
    if [ $iv1 -eq $iv2 ]; then continue; fi
    iv=$[10*$iv1 + $iv2]
    echo $iv1 $iv2 $iv
    echo $dir $modulation 2 $iv $phiR 1 >> $argfile
  done
done


# 3D
for ((iv1=0; iv1<$nIV; iv1++)); do
  if [ $allCombos -eq 1 ]; then start1=$[$iv1+1];
  else start1=0; fi
  for ((iv2=$start1; iv2<$nIV; iv2++)); do
    if [ $allCombos -eq 1 ]; then start2=$[$iv2+1];
    else start2=0; fi
    for ((iv3=$start2; iv3<$nIV; iv3++)); do
      if [ $iv1 -eq $iv2 ]; then continue; fi
      if [ $iv1 -eq $iv3 ]; then continue; fi
      if [ $iv2 -eq $iv3 ]; then continue; fi
      iv=$[100*$iv1 + 10*$iv2 + $iv3]
      echo $iv1 $iv2 $iv3 $iv
      echo $dir $modulation 3 $iv $phiR 1 >> $argfile
    done
  done
done


cat $argfile


# serial loop
while read args; do
  $exe $args
done < $argfile


mv -v spin*.root $outdir
mv -v *Canv_*.png $outdir
