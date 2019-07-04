#!/bin/bash

f=../mc/out_clasdispr.00.e10.600.emn0.75tmn.09.xs80.53nb.dis.0000.nrad.dat.evio.hipo
f=../skim/skim4_5036.hipo

make && valgrind \
  --leak-check=full \
  --log-file="outz" \
  --num-callers=40 \
  --suppressions=$ROOTSYS/etc/valgrind-root.supp \
  analysis.exe $f
  #--track-origins=yes \
