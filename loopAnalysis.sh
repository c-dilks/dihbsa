#!/bin/bash

if [ $# -ne 1 ]; then
  echo "usage: $0 [directory of hipo files]"
  exit
fi
indir=$1

job="batchAnalysis.bat";
log="logfiles"
mkdir -p $log

mkdir -p outroot/
rm -v spinroot/*.root
rm -v logfiles/*


echo "generate batchfile: $job"
echo "Executable = analysis.exe" > $job
echo "Universe = vanilla" >> $job
echo "notification = never" >> $job
echo "getenv = True" >> $job
echo "" >> $job


for file in ${indir}/*.hipo; do
  suffix=`echo $file|sed 's/^.*\///g'|sed 's/\.hipo$//g'`
  echo "prepare to analyze $file"
  echo "Arguments = $file" >> $job
  echo "Log    = ${log}/asym.${suffix}.log" >> $job
  echo "Output = ${log}/asym.${suffix}.out" >> $job
  echo "Error  = ${log}/asym.${suffix}.err" >> $job
  echo "Queue" >> $job
  echo "" >> $job
done

echo "submitting $job..."
condor_submit $job
