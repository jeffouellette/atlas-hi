# 2016 8TeV pPb main physics stream
./lib/runJetInsituCorrectionCheck 313063 0.03 true &
./lib/runJetInsituCorrectionCheck 313067 1.24 true &
./lib/runJetInsituCorrectionCheck 313100 9.66 true &
./lib/runJetInsituCorrectionCheck 313107 11.92 true & 
./lib/runJetInsituCorrectionCheck 313136 10.4 true & 
./lib/runJetInsituCorrectionCheck 313187 3.67 true &

wait

./lib/runJetInsituCorrectionCheck 313259 5.12 true &
./lib/runJetInsituCorrectionCheck 313285 4.74 true &
./lib/runJetInsituCorrectionCheck 313295 10.69 true &
./lib/runJetInsituCorrectionCheck 313333 4.13 true &
./lib/runJetInsituCorrectionCheck 313435 0.39 true &
./lib/runJetInsituCorrectionCheck 313572 0.01 false &

wait

./lib/runJetInsituCorrectionCheck 313574 1.33 false &
./lib/runJetInsituCorrectionCheck 313575 7.54 false &
./lib/runJetInsituCorrectionCheck 313603 8.69 false & 
./lib/runJetInsituCorrectionCheck 313629 6.86 false &
./lib/runJetInsituCorrectionCheck 313630 7.90 false & 
./lib/runJetInsituCorrectionCheck 313688 7.96 false &

wait

./lib/runJetInsituCorrectionCheck 313695 4.53 false & 
./lib/runJetInsituCorrectionCheck 313833 5.11 false &
./lib/runJetInsituCorrectionCheck 313878 2.16 false & 
./lib/runJetInsituCorrectionCheck 313929 0.63 false & 
./lib/runJetInsituCorrectionCheck 313935 10.96 false &
./lib/runJetInsituCorrectionCheck 313984 2.40 false &

wait

./lib/runJetInsituCorrectionCheck 314014 7.36 false & 
./lib/runJetInsituCorrectionCheck 314077 10.19 false &
./lib/runJetInsituCorrectionCheck 314105 6.50 false &
./lib/runJetInsituCorrectionCheck 314112 10.49 false &
./lib/runJetInsituCorrectionCheck 314157 9.83 false & 
./lib/runJetInsituCorrectionCheck 314170 4.92 false & 

wait

./lib/runJetInsituCorrectionCheckHist
