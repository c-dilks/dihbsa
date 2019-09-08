#!/bin/bash
# study asymmetry vs. run, or set of runs

# NOTE: must be executed with same ARGUMENTS you used when calling loopAsym.sh (except
# for the directory of ROOT files)

# loops through spinroot directory, linking $num runs at a time, computes their
# asymmetry, and renames spinFinal.root accordingly


##############
# $num is the number of runs to use to calculate each asymmetry; set this to 1
# if you want to calculate an asymmetry for every run
let num=4
##############



# move all spinroot files to a temp directory; the empty spinroot directory will 
# contain symlinks to root files in the temp directory
mkdir -p spinroot.tmp
mv spinroot{/*,.tmp/}

# initialize counters
let cnt=0
let tot=0
fname=""
nfiles=`ls spinroot.tmp|wc -l`


# MAIN LOOP
for file in spinroot.tmp/*.root; do

  # link a single file
  cd spinroot
  ln -s ../$file
  cd ..

  # get the run number
  runnum=`echo $file | sed 's/^.*spin\.skim._//g' | sed 's/\.hipo\.root$//g'`
  if [ $cnt -eq 0 ]; then fname="${fname}.${runnum}"; fi

  let cnt++
  let tot++

  # if we linked $num files, or if we've reached the end of the list of root files,
  # calculate the asymmetry and reset counters
  if [ $cnt -eq $num -o $tot -eq $nfiles ]; then

    echo "computing asymmetry for files:"
    ls -l --color=auto spinroot/

    # define file name which spinFinal.root will be renamed to
    fname="${fname}-${runnum}"

    # calcuate asymmetry
    asym.exe $* -c3
    sleep 3
    asym.exe $* -c4
    sleep 3
    mv spinFinal{,${fname}}.root

    # reset file names and clear out symlinks
    let cnt=0
    fname=""
    rm spinroot/*.root
    echo ""
  fi
done

# move root files from temp directory back to spinroot directory
mv spinroot{.tmp/*,/}
