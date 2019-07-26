#!/bin/bash

root -b -q PrintEnumerators.C > /dev/null
mkdir -p spinout

modulation=0
if [ $# -eq 0 ]; then
  echo "USAGE: $0 [outroot dir] [whichModulation]"
  exit
fi
if [ $# -ge 1 ]; then dir=$1; fi
if [ $# -ge 2 ]; then modulation=$2; fi


################
##### control booleans
allPermutations=0
dim1=1
dim2=0
dim3=0
################
################


argfile="args.dat"
> $argfile

# total number of IVs
nIV=4

# which phiR definition
phiR=3

# pairs list loop
while read pairTypeInfo; do
  echo "--- $pairTypeInfo"
  pairType=`echo $pairTypeInfo | awk '{print $1}'`
  pairName=`echo $pairTypeInfo | awk '{print $2}'`

# 1D
  if [ $dim1 -eq 1 ]; then
    for ((iv=0; iv<$nIV; iv++)); do
      echo $iv
      echo $dir $pairType $modulation 1 $iv $phiR 1 >> $argfile
    done
  fi


# 2D
  if [ $dim2 -eq 1 ]; then
    for ((iv1=0; iv1<$nIV; iv1++)); do
      if [ $allPermutations -eq 0 ]; then start=$[$iv1+1];
      else start=0; fi
      for ((iv2=$start; iv2<$nIV; iv2++)); do
        if [ $iv1 -eq $iv2 ]; then continue; fi
        iv=$[10*$iv1 + $iv2]
        echo $iv1 $iv2 $iv
        echo $dir $pairType $modulation 2 $iv $phiR 1 >> $argfile
      done
    done
  fi


# 3D
  if [ $dim3 -eq 1 ]; then
    for ((iv1=0; iv1<$nIV; iv1++)); do
      if [ $allPermutations -eq 0 ]; then start1=$[$iv1+1];
      else start1=0; fi
      for ((iv2=$start1; iv2<$nIV; iv2++)); do
        if [ $allPermutations -eq 0 ]; then start2=$[$iv2+1];
        else start2=0; fi
        for ((iv3=$start2; iv3<$nIV; iv3++)); do
          if [ $iv1 -eq $iv2 ]; then continue; fi
          if [ $iv1 -eq $iv3 ]; then continue; fi
          if [ $iv2 -eq $iv3 ]; then continue; fi
          iv=$[100*$iv1 + 10*$iv2 + $iv3]
          echo $iv1 $iv2 $iv3 $iv
          echo $dir $pairType $modulation 3 $iv $phiR 1 >> $argfile
        done
      done
    done
  fi

done < pairs.list


job="batchAsym.bat";
log="logfiles"
mkdir -p $log


echo "generate batchfile: $job"
echo "Executable = asym.exe" > $job
echo "Universe = vanilla" >> $job
echo "notification = never" >> $job
echo "getenv = True" >> $job
echo "" >> $job

while read args; do
  suffix=$(echo $args | sed 's/ /_/g')
  echo "Arguments = $args" >> $job
  echo "Log    = ${log}/asym.${suffix}.log" >> $job
  echo "Output = ${log}/asym.${suffix}.out" >> $job
  echo "Error  = ${log}/asym.${suffix}.err" >> $job
  echo "Queue" >> $job
  echo "" >> $job
done < $argfile

echo "submitting $job..."
condor_submit $job
sleep 1
condor_q
