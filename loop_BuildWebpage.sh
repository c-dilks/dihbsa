#!/bin/bash

htmlfile="spinIndex.html"
> $htmlfile
function a { 
  echo "$1" >> $htmlfile
}

a "<html><head>"
a "<title>spin asymmetries</title>"
a "</head><body>"

for dir in `ls -d spinout*/`; do
  ./BuildWebpage.sh $dir
  html=$(ls -t *.html|head -n1)
  name=$(echo $html | sed 's/^spin\.//' | sed 's/\.html$//')
  a "<a href=\"${html}\">${name}</a><br />"
done

a "</body></html>"
