#!/bin/bash

for fileS in spinout/kindepCanv_sinPhiHR_scale_*; do 
  fileW=$(echo $fileS|sed 's/_scale/_weight/g')
  fileN=$(echo $fileS|sed 's/_scale//g')
  modN=$(echo $fileN|sed 's/kindep/asymMod/')
  modS=$(echo $fileS|sed 's/kindep/asymMod/')
  modW=$(echo $fileW|sed 's/kindep/asymMod/')
  if [ -f $modN ]; then
    sxiv $fileN $fileW $fileS $modN $modW $modS
  else
    sxiv $fileN $fileW $fileS
  fi
done
