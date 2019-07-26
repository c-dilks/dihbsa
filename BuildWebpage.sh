#!/bin/bash

pairStr="piPlus_piMinus"
if [ $# -ne 2 ]; then 
  echo "USAGE: $0 [pairStr] [spinout_dir] "
  echo "... prefer to call this by loop_BuildWebpage.sh"
  exit
fi
pairStr=$1
dir=$2

html="${dir}/web_${pairStr}.html"

function a { 
  echo "$1" >> $html
}

function ai {
  a "<a href=\"${1}_${2}.${pairStr}.png\">"
  a "<img class=\"${3}\" src=\"${1}_${2}.${pairStr}.png\" />"
  a "</a>"
}


> $html
> suffixList1
> suffixList2
> suffixList3


for f in ${dir}/kindepCanv*.${pairStr}.png; do
  suf=$(echo $f | sed 's/^.*\/kindepCanv_//g' | sed 's/\..*\.png//g')
  echo $suf

  if [ -z `echo $suf|grep bins` ]; then 
    # 1D
    vars=$(echo $suf| sed 's/^.*_//g')
    echo $suf >> suffixList1
  else
    wordlist=$(echo $suf | sed 's/_/ /g')
    wordcnt=$(echo $wordlist | wc -w)

    let cnt=$wordcnt
    w1=$(echo $wordlist | cut -d" " -f$cnt)
    let cnt--
    w2=$(echo $wordlist | cut -d" " -f$cnt)
    let cnt--
    w3=$(echo $wordlist | cut -d" " -f$cnt)
    let cnt--
    w4=$(echo $wordlist | cut -d" " -f$cnt)

    if [ "$w2" == "bins" ]; then
      echo $suf >> suffixList2
    elif [ "$w3" == "bins" ]; then
      echo $suf >> suffixList3
    fi

  fi
done

echo 1D:
cat suffixList1
echo 2D:
cat suffixList2
echo 3D:
cat suffixList3


a "<html><head>"
a "<title>${pairStr} asymmetries</title>"
a "<style>"
a ".h { height: 400px; }"
a ".w { width: 800px; }"
a "</style>"
a "</head><body>"
a "<h1>${pairStr} asymmetries</h1><hr />"


a "<br /><h1>1D binning:</h1><hr />"
while read s; do
  for ((amp=0; amp<2; amp++)); do ai RF_A${amp}_kindepCanv $s h; done
  a "<br />"
  ai kindepCanv $s h
  ai asymModCanv $s h
  a "<br />"
  ai chindfCanv $s h
  ai modDistCanv $s h
  a "<br />"
  ai rellumCanv $s h
  a "<hr />"
done < suffixList1


a "<br /><br /><h1>2D binning:</h1><hr />"
while read s; do
  for ((amp=0; amp<2; amp++)); do ai RF_A${amp}_kindepCanv $s h; done
  a "<br />"
  ai kindepCanv $s h
  ai asymModCanv $s h
  a "<br />"
  ai modDistCanv $s h
  a "<br />"
  ai chindfCanv $s h
  a "<br />"
  ai rellumCanv $s h
  a "<hr />"
done < suffixList2


a "<br /><br /><h1>3D binning:</h1><hr />"
while read s; do
  ai kindepCanv $s w
  ai chindfCanv $s w
  a "<br />"
  ai rellumCanv $s w
  a "<hr />"
done < suffixList3


a "</body></html>"


rm suffixList1
rm suffixList2
rm suffixList3

echo "---"
echo ${html} built
