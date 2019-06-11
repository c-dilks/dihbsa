#!/bin/bash

root -b -q PrintEnumerators.C

dir="outroot.dnp2018.old"

while read line; do
  pairtype=`echo $line | awk '{print $1}'`
  pairname=`echo $line | awk '{print $2}'`
  loop_asym.sh $dir "$pairtype"
  sleep 3
  mv -v spinout{,.${pairname}}
  sleep 3
done < pairs.list
