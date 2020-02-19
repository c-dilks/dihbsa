#!/bin/bash
for p in chris tim; do
  for f in spinroot.${p}Grids/asym*.root; do
    root -b -q drawAsymmetries.C'("'${f}'")'
  done
done
