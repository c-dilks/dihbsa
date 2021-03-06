#!/bin/bash

dir="spinout"
if [ $# -gt 0 ]; then dir="$1"; fi

for file in ${dir}/*.root; do
  root -b -q drawAsymPlots.C'("'$file'")'
done
exit # imagemagick is currently broken?

# imagemagick
for file in ${dir}/*multi*.pdf; do
  echo "convert $file to png..."
  convert -trim -density 300 ${file} `echo $file|sed 's/pdf$/png/'`
done
