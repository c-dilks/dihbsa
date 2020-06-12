#!/bin/bash
#calls Projector.C for various cases

# target info
tarTitle=('NH_{3}' 'ND_{3}' '^{3}He')
tarName=( 'NH3'    'ND3'    'He3')
tarPol=(  0.85     0.35     0.50) # polarization
tarDil=(  0.20     0.285    0.27) # dilution
modelH=(  'xdepP'  'xdepN'  'xdepN') # model for TSA hL(x) projection
modelG=(  'specP'  'specN'  'specN') # model for TSA G1perp(Mh) projection
modelD=(  'zero'   'zero'   'zero') # model for DSA projection

# beam polarization (for DSA only)
beamPol=0.8

# subroutine to rename output file and build tables
function organize {
# args: 1=rootName 2=tarName 3=days 4={H,G,DSA}
  mv -v ${1}.proj.pdf projection${4}_${2}_${3}days.pdf
  if [ $3 -eq 60 ]; then
    cp binranges.tabletex projection${4}_${2}.tabletex
  fi
  paste projection${4}_${2}.tabletex ${1}.tabletex > tempo
  mv tempo projection${4}_${2}.tabletex
}


# loop over number of days
for days in 60 90 120; do

  # loop over targets
  for i in {0..2}; do


# projection for hL(x)
    root -b -q Projector.C\
'("spinroot.proj.vsX/asym_4.root",'\
'"'${modelH[$i]}'",'\
'"'${tarTitle[$i]}'",'\
${days}','\
1','\
${tarPol[$i]}','\
${tarDil[$i]}')'
  organize kindepMA_A2_X ${tarName[$i]} ${days} H

# projection for G1perp
    root -b -q Projector.C\
'("spinroot.proj.vsM/asym_4.root","'\
${modelG[$i]}'","'\
${tarTitle[$i]}'",'\
${days}','\
1','\
${tarPol[$i]}','\
${tarDil[$i]}')'
  organize kindepMA_A1_M ${tarName[$i]} ${days} G

# projection for DSA
    root -b -q Projector.C\
'("spinroot.proj.vsZ/asym_1000.root","'\
${modelD[$i]}'","'\
${tarTitle[$i]}'",'\
${days}','\
${beamPol}','\
${tarPol[$i]}','\
${tarDil[$i]}')'
  organize kindepMA_A1_Z ${tarName[$i]} ${days} DSA

  done
done

rm kindep*.tabletex
rm binranges*.tabletex


# add tabletex headers
function addHeaders {
# args: 1={M_h,x,z} 2={UL,LL} 3={H,G,DSA}
  sigmaStr=""
  sigmaStr='60 days $\sigma$ & 90 days $\sigma$ & 120 days $\sigma$ \\\hline'
  for file in projection${3}*.tabletex; do
    
    fn=$(echo $file|sed 's/\./_/g')

    echo '\begin{tabular}{|c|c|c|c|c|c|}' >> tempo
    echo '\hline' >> tempo
    echo '$'${1}'$ range & $\langle '${1}' \rangle$ & $A_{'${2}'}$ & '$sigmaStr >> tempo

    while read line; do echo ${line}' \\\hline' >> tempo; done < $file

    echo '\end{tabular}' >> tempo
    echo '\caption{'$fn'}' >> tempo
    echo '\label{tab:'$fn'}' >> tempo
    echo '' >> tempo

    mv tempo $file
  done
}
addHeaders x UL H
addHeaders M_h UL G
addHeaders z LL DSA
