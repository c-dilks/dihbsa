#!/bin/bash
# loops MCmatch.C over all files in outroot.MC.{gen,rec}
# outputs histograms which indicate how well MCgen and MCrec events are matched

touch match.root
rm match.root

let limiter=0
for file in outroot.MC.gen/*.root; do
  echo "analyze $file"
  root -b -q MCmatch.C'("'$(echo $file | sed 's/^.*\///g')'")'
  let limiter++
  #if [ $limiter -gt 30 ]; then break; fi
done

hadd tmp.root match*.root
rm match.*.root
mv {tmp,match}.root

root -l drawMatchFraction.C

#cat > tmp.C << 'EOF'
#void tmp() {
  #TFile * f = new TFile("match.root","READ");
  #TTree * mtr = (TTree*) f->Get("mtr");
  #//mtr->Draw("diff_hadE[0]:diff_hadE[1]>>a(100,0,0.1,100,0,0.1)","","colz");
#}
#EOF

#root -l tmp.C
#rm tmp.C
