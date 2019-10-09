#!/bin/bash

whichformu="three"

pairtype="0x34"
mod=2

mkdir -p forPrelim
root -b -q PrintEnumerators.C

while read pt; do pairtype=`echo $pt|awk '{print $1}'`
#for mod in 0 2; do
for iv in 1 2 3 32; do
  dest="forPrelim/m${mod}.${pairtype}.i${iv}.${whichformu}"
  loopAsym.sh outroot.fall18 -p $pairtype -m$mod -i$iv \
    2>&1 | tee ${dest}.log
  sleep 3
  echo $dest
  mv -v spinFinal.root ${dest}.root
  sleep 1
done
#done
done < pairs.list
