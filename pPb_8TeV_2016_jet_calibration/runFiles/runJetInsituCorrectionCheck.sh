# 2016 8TeV pPb main physics stream
./bin/runJetInsituCorrectionCheck 313063 0.03 true &
./bin/runJetInsituCorrectionCheck 313067 1.24 true &
./bin/runJetInsituCorrectionCheck 313100 9.66 true &
./bin/runJetInsituCorrectionCheck 313107 11.92 true & 
./bin/runJetInsituCorrectionCheck 313136 10.4 true & 
./bin/runJetInsituCorrectionCheck 313187 3.67 true &

wait

./bin/runJetInsituCorrectionCheck 313259 5.12 true &
./bin/runJetInsituCorrectionCheck 313285 4.74 true &
./bin/runJetInsituCorrectionCheck 313295 10.69 true &
./bin/runJetInsituCorrectionCheck 313333 4.13 true &
./bin/runJetInsituCorrectionCheck 313435 0.39 true &
./bin/runJetInsituCorrectionCheck 313572 0.01 false &

wait

./bin/runJetInsituCorrectionCheck 313574 1.33 false &
./bin/runJetInsituCorrectionCheck 313575 7.54 false &
./bin/runJetInsituCorrectionCheck 313603 8.69 false & 
./bin/runJetInsituCorrectionCheck 313629 6.86 false &
./bin/runJetInsituCorrectionCheck 313630 7.90 false & 
./bin/runJetInsituCorrectionCheck 313688 7.96 false &

wait

./bin/runJetInsituCorrectionCheck 313695 4.53 false & 
./bin/runJetInsituCorrectionCheck 313833 5.11 false &
./bin/runJetInsituCorrectionCheck 313878 2.16 false & 
./bin/runJetInsituCorrectionCheck 313929 0.63 false & 
./bin/runJetInsituCorrectionCheck 313935 10.96 false &
./bin/runJetInsituCorrectionCheck 313984 2.40 false &

wait

./bin/runJetInsituCorrectionCheck 314014 7.36 false & 
./bin/runJetInsituCorrectionCheck 314077 10.19 false &
./bin/runJetInsituCorrectionCheck 314105 6.50 false &
./bin/runJetInsituCorrectionCheck 314112 10.49 false &
./bin/runJetInsituCorrectionCheck 314157 9.83 false & 
./bin/runJetInsituCorrectionCheck 314170 4.92 false & 

wait

./bin/runJetInsituCorrectionCheckHist
