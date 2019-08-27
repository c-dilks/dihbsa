#!/bin/bash

root -b -q PrintEnumerators.C

dir="outroot.dnp2018"
modulation=5
if [ $# -ge 1 ]; then dir=$1; fi
#if [ $# -ge 2 ]; then modulation=$2; fi

#for modulation in 0 3; do
  while read line; do
    pairtype=`echo $line | awk '{print $1}'`
    pairname=`echo $line | awk '{print $2}'`
    loop_asym.sh $dir "$pairtype" $modulation
    sleep 3
    mv -v spinout{,.${pairname}}
    sleep 3
  done < pairs.list

  mkdir -p spinweb${modulation}
  mv spinout.* spinweb${modulation}/
#done
