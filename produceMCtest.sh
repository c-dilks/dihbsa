#!/bin/bash
for h in {0..11}; do
  for data in gen rec; do
    [[ $h -lt 10 ]] && hh="0$h" || hh="$h"
    dest="forMC/${data}.${hh}"
    loopAsym.sh outroot.MC.${data} -m2 -h$h 2>&1 | tee ${dest}.log
    mv -v {spinFinal,$dest}.root
  done
done
