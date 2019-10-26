#!/bin/bash

dir="forMC"

for file in ${dir}/*.root; do
  root -b -q drawAsymPlots.C'("'$file'")'
done

# imagemagick
for file in ${dir}/*multi*.pdf; do
  echo "convert $file to png..."
  convert -trim -density 300 ${file} `echo $file|sed 's/pdf$/png/'`
done
