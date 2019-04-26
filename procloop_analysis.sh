#!/bin/bash

exe="analysis.exe"

# default arguments
hipodir=../hipo  # directory of hipo files to analyse
njobsMax=4  # maximum number of parallel jobs to run
            # (set this to zero to use all accessible threads)

if [ $# -ge 1 ]; then hipodir=$1; fi
if [ $# -ge 2 ]; then njobsMax=$2; fi


# determine number of accessible threads and
# based on that, determine number of parallel jobs to run
nthreads=$(grep -c processor /proc/cpuinfo)
if [ $njobsMax -eq 0 ]; then
  njobs=$nthreads
elif [ $njobsMax -gt $nthreads ]; then
  njobs=$nthreads
else
  njobs=$njobsMax
fi

# print out config
echo "[+] reading HIPO files from $hipodir"
echo "[+] found $nthreads accessible CPU threads"
echo "[+] max number of parallel jobs to run: $njobs"


# make subdirs
logdir="logfiles"
outdir="outroot"
mkdir -p $logdir
mkdir -p $outdir


# filter for hipo files
#hipofilter="32.evio.8"
hipofilter="evio"

# get list of hipo files
ls ${hipodir}/*.hipo | grep $hipofilter > hipoz
filecnt=$(cat hipoz | wc -l)
echo "[+] number of hipo files to process: $filecnt"


# execute analysis code in parallel
cnt=1
echo ""
echo "[+] begin job submission"
while read hipofile; do
  log=$(echo $hipofile | sed 's/^.*\///g' | sed 's/hipo$/log/g')
  if [ $cnt -le $njobs ]; then 
    ./$exe $hipofile >& $logdir/Phi$log &
    echo " processing $hipofile"
    let cnt++
  else
    wait
    cnt=1
  fi
done < hipoz
rm hipoz


# echo final output
echo "[+] done submitting all jobs"
echo ""
echo "[+] list of jobs still running in background:"
#ps -u `whoami` -f | grep $exe
pgrep -u `whoami` -a $exe
echo ""
echo "[+] use the following command to check on them:"
#echo "ps -u `whoami` -f | grep $exe"
echo "pgrep -u `whoami` -a $exe"

