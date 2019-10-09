#!/bin/bash
#for file in forMC/*.root; do
for file in forPrelim/*.root; do
  root -b -q drawPrelim.C'("'$file'")'
done
