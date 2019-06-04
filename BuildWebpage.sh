#!/bin/bash

dir="spinout.pairP0"
if [ $# -ge 1 ]; then dir="$1"; fi

name=$(echo $dir | sed 's/spinout\.//g' | sed 's/\/$//g')
html=spin.${name}.html

function a { 
  echo "$1" >> $html
}

function ai {
  a "<a href=\"${dir}/${1}_${2}.png\">"
  a "<img class=\"${3}\" src=\"${dir}/${1}_${2}.png\" />"
  a "</a>"
}


> $html
> suffixList1
> suffixList2
> suffixList3

for f in ${dir}/kindepCanv*.png; do
  suf=$(echo $f | sed 's/^.*\/kindepCanv_//g' | sed 's/\.png//g')

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
a "<title>$name</title>"
a "<style>"
a ".h { height: 400px; }"
a ".w { width: 800px; }"
a "</style>"
a "</head><body>"
a "<h1>${name}</h1><hr />"


a "<br /><h1>1D:</h1><hr />"
while read s; do
  ai kindepCanv $s h
  ai asymModCanv $s h
  a "<br />"
  ai chindfCanv $s h
  ai rellumCanv $s h
  a "<hr />"
done < suffixList1


a "<br /><br /><h1>2D:</h1><hr />"
while read s; do
  ai kindepCanv $s h
  ai asymModCanv $s h
  a "<br />"
  ai chindfCanv $s h
  a "<br />"
  ai rellumCanv $s h
  a "<hr />"
done < suffixList2


a "<br /><br /><h1>3D:</h1><hr />"
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
