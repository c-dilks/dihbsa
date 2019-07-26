#!/bin/bash

if [ $# -ne 1 ]; then 
  echo "USAGE: $0 [spinout_dir]"
  exit
fi
dir=$1

modulation=$(echo $dir | sed 's/^spinout_//g')
htmlfile="spinIndex_${modulation}.html"
> $htmlfile

function a { 
  echo "$1" >> $htmlfile
}

function mkHadStr {
  if [ "$1" == "$2" ]; then hadStr="${1}1_${2}2"
  else hadStr="${1}_${2}"; fi
}


a "<html><head>"
a "<title>spin asymmetries -- ${modulation}</title>"
a "</head><body>"

root -b -q PrintEnumerators.C'(true)'

hadList="piPlus piMinus diphoton KPlus KMinus"

a "<h1>${modulation} asymmetries</h1><hr />"
a "<table style='width:80%; height:80%'>"
a "<th>"
for had in $hadList; do a "<td><b>${had}</b></td>"; done
a "</th>"
for had1 in $hadList; do
  a "<tr>"
  a "<td><b>$had1</b></td>"
  for had2 in $hadList; do
    a "<td>"

    # build hadron pair name; if it's not in pairs.list, then the hadron ordering
    # is switched
    mkHadStr $had1 $had2
    if [ -z "`grep $hadStr pairs.list`" ]; then mkHadStr $had2 $had1; fi

    echo "$had1 $had2 -- $hadStr"
    a "<a href=\"${dir}/web_${hadStr}.html\">${had1}&nbsp;&nbsp;&nbsp;${had2}</a>"
    ./BuildWebpage.sh $hadStr $dir

    a "</td>"
  done
  a "</tr>"
done
a "</table>"
a "</body></html>"
