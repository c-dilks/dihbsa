#!/bin/bash

#hipoDir="/w/hallb-scifs17exp/clas12/rg-a/data/cooked_5p7p8_fullMap_alignprod"
hipoDir="/home/dilks/j/hipo"
workDir=$(pwd -P)

ls $hipoDir | sed 's/\.hipo$//' > hipoFiles.list


outxml="test.xml"
> $outxml


function app {
  echo "$1" >> $outxml
}

app "<Request>"
app "<Project name=\"clas12\" />"
app "<Track name=\"debug\" />"
app "<Name name=\"dihbsaAnalysis\" />"
app ""

app "<Variable name=\"hipoDir\" value=\"file:${hipoDir}\" />"
app "<Variable name=\"workDir\" value=\"file:${workDir}\" />"
app ""

app "<List name=\"hipoFiles\">"
while read hipo; do app "$hipo"; done < hipoFiles.list
app "</List>"
app ""

app "<ForEach list=\"hipoFiles\">"
app "<Job>"
app "  <TimeLimit unit=\"minutes\" time=\"20\" />"
app "  <DiskSpace space=\"2\" unit=\"GB\" />"
app "  <Memory space=\"2\" unit=\"GB\" />"
app ""

app "  <Input src=\"\${hipoDir}/\${hipoFiles}.hipo\" dest=\"infile.hipo\" />"
app ""

app "  <Command><![CDATA["
app "    analysis.exe infile.hipo"
app "  ]]></Command>"
app ""

app "  <Output dest=\"\${workDir}/outroot/\${hipoFiles}.root\" src=\"outroot.root\" />"
app "  <Stdout dest=\"\${workDir}/logfiles/\${hipoFiles}.out\" />"
app "  <Stderr dest=\"\${workDir}/logfiles/\${hipoFiles}.err\" />"
app ""

app "</Job>"
app "</ForEach>"
app ""
app "</Request>"


echo "----------------"
cat $outxml
