#!/bin/bash
# build config file for EIC

if [ $# -ne 2 ]; then
  echo "USAGE: $0 EbeamEn PbeamEn"
  exit
fi
> config.cf
function app { echo "$1 = $2" >> config.cf; }

app Experiment eic
app EbeamEn $1
app PbeamEn $2

cat config.cf
