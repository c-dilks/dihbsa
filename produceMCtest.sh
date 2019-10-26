#!/bin/bash
for h in {0..10}; do
  for data in gen rec; do
    dest="forMC/${data}.${h}"
    loopAsym.sh outroot.MC.${data} -m2 -h$h 2>&1 | tee ${dest}.log
    mv -v {spinFinal,$dest}.root
  done
done
