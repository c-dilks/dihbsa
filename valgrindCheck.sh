#!/bin/bash

make && valgrind \
  --leak-check=full \
  --log-file="outz" \
  --num-callers=40 \
  --suppressions=$ROOTSYS/etc/valgrind-root.supp \
  catSpinroot.exe
  #asym.exe -p 0x34 -m2 -i2 -c3
  #--track-origins=yes \
