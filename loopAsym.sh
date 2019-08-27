#!/bin/bash

if [ $# -eq 0 ]; then
  echo "usage: $0 [ARGUMENTS]"
  echo "the first argument must be the directory of ROOT files to analyze"
  echo "the remaining arguments can be any of the OPTIONS listed when you"
  echo "execute asym.exe without any arguments"
  echo "WARNING: do NOT use the flow control variable when executing this script!"
  exit
fi

args=$*
indir=`echo $args|awk '{print $1}'`
opts=`echo $args|sed 's/[^ ]* *//'`

job="batchAsym.bat";
log="logfiles"
mkdir -p $log

mkdir -p spinroot

echo "generate batchfile: $job"
echo "Executable = asym.exe" > $job
echo "Universe = vanilla" >> $job
echo "notification = never" >> $job
echo "getenv = True" >> $job
echo "" >> $job


for file in ${indir}/*.root; do
  suffix=`echo $file|sed 's/^.*\///g'|sed 's/\.root$//g'`
  echo "prepare to analyze $file"
  echo "Arguments = -f $file $opts -c2" >> $job
  echo "Log    = ${log}/asym.${suffix}.log" >> $job
  echo "Output = ${log}/asym.${suffix}.out" >> $job
  echo "Error  = ${log}/asym.${suffix}.err" >> $job
  echo "Queue" >> $job
  echo "" >> $job
done

#exit

echo "submitting $job..."
condor_submit $job
sleep 1
condor_q

./wait_for_condor
asym.exe -c3
