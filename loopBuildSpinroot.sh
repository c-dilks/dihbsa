#!/bin/bash

if [ $# -eq 0 ]; then
  echo "usage: $0 [ARGUMENTS]"
  echo "the first argument must be the directory of outroot files to analyze"
  echo "the remaining arguments can be any of the options listed when you"
  echo "execute buildSpinroot.exe without any arguments"
  exit
fi

args=$*
indir=`echo $args|awk '{print $1}'`
opts=`echo $args|sed 's/[^ ]* *//'`

job="batchBuildSpinroot.bat";
log="logfiles"
mkdir -p $log

mkdir -p spinroot
rm -v spinroot/*.root


echo "generate batchfile: $job"
echo "Executable = buildSpinroot.exe" > $job
echo "Universe = vanilla" >> $job
echo "notification = never" >> $job
echo "getenv = True" >> $job
echo "" >> $job


for file in ${indir}/*.root; do
  suffix=`echo $file|sed 's/^.*\///g'|sed 's/\.root$//g'`
  echo "prepare to analyze $file"
  echo "Arguments = -f $file $opts" >> $job
  echo "Log    = ${log}/buildSpinroot.${suffix}.log" >> $job
  echo "Output = ${log}/buildSpinroot.${suffix}.out" >> $job
  echo "Error  = ${log}/buildSpinroot.${suffix}.err" >> $job
  echo "Queue" >> $job
  echo "" >> $job
done

echo "submitting $job..."
condor_submit $job
