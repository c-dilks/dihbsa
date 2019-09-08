#!/bin/bash
for file in spinFinal.5*.root; do
  root -b -q DrawMultiGr.C'("'$file'")'
done
