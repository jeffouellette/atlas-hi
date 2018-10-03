# 2016 8TeV pPb main physics stream
./lib/runZmumuJets 313063 0.03 false true &
./lib/runZmumuJets 313067 1.24 false true &
./lib/runZmumuJets 313100 9.66 false true &
./lib/runZmumuJets 313107 11.92 false true & 
./lib/runZmumuJets 313136 10.4 false true & 
./lib/runZmumuJets 313187 3.67 false true &

wait

./lib/runZmumuJets 313259 5.12 false true &
./lib/runZmumuJets 313285 4.74 false true &
./lib/runZmumuJets 313295 10.69 false true &
./lib/runZmumuJets 313333 4.13 false true &
./lib/runZmumuJets 313435 0.39 false true &
./lib/runZmumuJets 313572 0.01 false false &

wait

./lib/runZmumuJets 313574 1.33 false false &
./lib/runZmumuJets 313575 7.54 false false &
./lib/runZmumuJets 313603 8.69 false false & 
./lib/runZmumuJets 313629 6.86 false false &
./lib/runZmumuJets 313630 7.90 false false & 
./lib/runZmumuJets 313688 7.96 false false &

wait

./lib/runZmumuJets 313695 4.53 false false & 
./lib/runZmumuJets 313833 5.11 false false &
./lib/runZmumuJets 313878 2.16 false false & 
./lib/runZmumuJets 313929 0.63 false false & 
./lib/runZmumuJets 313935 10.96 false false &
./lib/runZmumuJets 313984 2.40 false false &

wait

./lib/runZmumuJets 314014 7.36 false false & 
./lib/runZmumuJets 314077 10.19 false false &
./lib/runZmumuJets 314105 6.50 false false &
./lib/runZmumuJets 314112 10.49 false false &
./lib/runZmumuJets 314157 9.83 false false & 
./lib/runZmumuJets 314170 4.92 false false & 

wait

# Zmumu (data overlay) samples for period A
./lib/runZmumuJets 0 0 true true user.jeouelle.2.4.30hi.calibcheck.200.mc15_8TeV.361107.ZmumuJet.pPb_myOutput.root &

# ... and for period B
./lib/runZmumuJets 0 0 true false user.jeouelle.2.4.30hi.calibcheck.200.mc15_8TeV.361107.ZmumuJet.Pbp_myOutput.root &

wait

./lib/runZmumuJetsHist  
