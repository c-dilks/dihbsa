#!/bin/bash

root -b -q PrintEnumerators.C

dir="outroot.dnp2018"

mkdir -p diagset

while read line; do
  pairtype=`echo $line | awk '{print $1}'`
  pairname=`echo $line | awk '{print $2}'`
  diagnostics.exe $dir $pairtype
  sleep 3
  mv -v plots.root diagset/plots.${pairname}.root
  sleep 3
done < pairs.list
