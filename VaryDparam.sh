#!/bin/bash

# fraction dihadrons from Vector Mesons (evaluated in MC)
# - this is used as D_{1,OO}^{p} / D_{1,OO} in the positivity bounds
VMfrac=0.12
echo "VMfrac = ${VMfrac}"

# positivity bounds
lb=$(echo "-3.0/2.0*${VMfrac}" | bc -l)
ub=$(echo "3.0*${VMfrac}" | bc -l)
echo "D_{1,LL} / D_{1,OO} will be varied between ${lb} and ${ub}"

# determine sequence of Dparam values
nsteps=5
stepsize=$(echo "(${ub}-${lb})/${nsteps}" | bc -l)
#seq $lb $stepsize $ub

