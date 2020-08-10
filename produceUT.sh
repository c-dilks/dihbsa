asymFit.exe 921 0 | tee fit_test1.log
mv spinroot/asym_{921,test1}.root

asymFit.exe 801 0 | tee fit_test2.log
mv spinroot/asym_{801,test2}.root

asymFit.exe 803 0 | tee fit_test3.log
mv spinroot/asym_{803,test3}.root

mv fit_test*.log spinroot/
