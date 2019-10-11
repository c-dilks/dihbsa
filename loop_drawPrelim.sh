#!/bin/bash
for file in forMC/*.root; do
for file in forPrelim/*.root; do
  root -b -q drawPrelim.C'("'$file'")'
done

# imagemagick
for file in forPrelim/*multi*.pdf; do
  echo "convert $file to png..."
  convert -trim -density 300 ${file} `echo $file|sed 's/pdf$/png/'`
done
