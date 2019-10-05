#!/bin/bash

pairtype="0x34"
mod=2
whichformu="one"

mkdir -p forPrelim

for iv in {1..3}; do
  #loopAsym.sh outroot.fall18 -p $pairtype -m$mod -i$iv
  #sleep 3
  dest="forPrelim/m${mod}.${pairtype}.i${iv}.${whichformu}.root"
  echo dest
  #mv -v spinFinal.root $dest
  #sleep 1
done
