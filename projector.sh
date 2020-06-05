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

# loop over number of days
for days in 60 90 120; do

  # loop over targets
  for i in {0..2}; do

    suffix="${tarName[$i]}_${days}days"

# projection for hL(x)
    root -b -q Projector.C\
'("spinroot.proj.vsX/asym_4.root",'\
'"'${modelH[$i]}'",'\
'"'${tarTitle[$i]}'",'\
${days}','\
1','\
${tarPol[$i]}','\
${tarDil[$i]}')'
    mv -v kindepMA_A2_X.proj.pdf projectionH_${suffix}.pdf

# projection for G1perp
    root -b -q Projector.C\
'("spinroot.proj.vsM/asym_4.root","'\
${modelG[$i]}'","'\
${tarTitle[$i]}'",'\
${days}','\
1','\
${tarPol[$i]}','\
${tarDil[$i]}')'
    mv -v kindepMA_A1_M.proj.pdf projectionG_${suffix}.pdf

# projection for DSA
    root -b -q Projector.C\
'("spinroot.proj.vsZ/asym_1000.root","'\
${modelD[$i]}'","'\
${tarTitle[$i]}'",'\
${days}','\
${beamPol}','\
${tarPol[$i]}','\
${tarDil[$i]}')'
    mv -v kindepMA_A1_Z.proj.pdf projectionDSA_${suffix}.pdf

  done
done


# projection for DSA partial waves
##root -b -q Projector.C\
##'("spinroot.proj.vsZ/asym_1001.root","zero",30,0.85,0.85,0.2)'
##mv kindepMA_A0_Z.proj.pdf projectionPW0.pdf
##mv kindepMA_A1_Z.proj.pdf projectionPW1.pdf
##mv kindepMA_A2_Z.proj.pdf projectionPW2.pdf
