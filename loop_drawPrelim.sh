#!/bin/bash
for file in forPrelim/*.root; do
  root -b -q drawPrelim.C'("'$file'")'
done
