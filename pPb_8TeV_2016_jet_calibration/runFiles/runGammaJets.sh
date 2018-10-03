# 2016 8TeV pPb main physics stream
./lib/runGammaJets 313063 0.03 false true &
./lib/runGammaJets 313067 1.24 false true &
./lib/runGammaJets 313100 9.66 false true &
./lib/runGammaJets 313107 11.92 false true & 
./lib/runGammaJets 313136 10.4 false true & 
./lib/runGammaJets 313187 3.67 false true &

wait

./lib/runGammaJets 313259 5.12 false true &
./lib/runGammaJets 313285 4.74 false true &
./lib/runGammaJets 313295 10.69 false true &
./lib/runGammaJets 313333 4.13 false true &
./lib/runGammaJets 313435 0.39 false true &
./lib/runGammaJets 313572 0.01 false false &

wait

./lib/runGammaJets 313574 1.33 false false &
./lib/runGammaJets 313575 7.54 false false &
./lib/runGammaJets 313603 8.69 false false & 
./lib/runGammaJets 313629 6.86 false false &
./lib/runGammaJets 313630 7.90 false false & 
./lib/runGammaJets 313688 7.96 false false &

wait

./lib/runGammaJets 313695 4.53 false false & 
./lib/runGammaJets 313833 5.11 false false &
./lib/runGammaJets 313878 2.16 false false & 
./lib/runGammaJets 313929 0.63 false false & 
./lib/runGammaJets 313935 10.96 false false &
./lib/runGammaJets 313984 2.40 false false &

wait

./lib/runGammaJets 314014 7.36 false false & 
./lib/runGammaJets 314077 10.19 false false &
./lib/runGammaJets 314105 6.50 false false &
./lib/runGammaJets 314112 10.49 false false &
./lib/runGammaJets 314157 9.83 false false & 
./lib/runGammaJets 314170 4.92 false false & 

wait

# gamma+jet pp sample with pPb data overlay
./lib/runGammaJets 1 0 true true user.jeouelle.2.4.30hi.calibcheck.200.mc15_8TeV.423100.Slice1.pPb_myOutput.root &
./lib/runGammaJets 2 0 true true user.jeouelle.2.4.30hi.calibcheck.200.mc15_8TeV.423101.Slice2.pPb_myOutput.root &
./lib/runGammaJets 3 0 true true user.jeouelle.2.4.30hi.calibcheck.200.mc15_8TeV.423102.Slice3.pPb_myOutput.root &
./lib/runGammaJets 4 0 true true user.jeouelle.2.4.30hi.calibcheck.200.mc15_8TeV.423103.Slice4.pPb_myOutput.root &
./lib/runGammaJets 5 0 true true user.jeouelle.2.4.30hi.calibcheck.200.mc15_8TeV.423104.Slice5.pPb_myOutput.root &
./lib/runGammaJets 6 0 true true user.jeouelle.2.4.30hi.calibcheck.200.mc15_8TeV.423105.Slice6.pPb_myOutput.root &

wait

./lib/runGammaJets 1 0 true false user.jeouelle.2.4.30hi.calibcheck.200.mc15_8TeV.423100.Slice1.Pbp_myOutput.root &
./lib/runGammaJets 2 0 true false user.jeouelle.2.4.30hi.calibcheck.200.mc15_8TeV.423101.Slice2.Pbp_myOutput.root &
./lib/runGammaJets 3 0 true false user.jeouelle.2.4.30hi.calibcheck.200.mc15_8TeV.423102.Slice3.Pbp_myOutput.root &
./lib/runGammaJets 4 0 true false user.jeouelle.2.4.30hi.calibcheck.200.mc15_8TeV.423103.Slice4.Pbp_myOutput.root &
./lib/runGammaJets 5 0 true false user.jeouelle.2.4.30hi.calibcheck.200.mc15_8TeV.423104.Slice5.Pbp_myOutput.root &
./lib/runGammaJets 6 0 true false user.jeouelle.2.4.30hi.calibcheck.200.mc15_8TeV.423105.Slice6.Pbp_myOutput.root &

wait

# alternative validation signal-only gamma+jet sample (period A only!)
./lib/runGammaJets 1 0 true true user.jeouelle.2.4.30hi.calibcheck.signalonly.310.mc15_8TeV.423100.Slice1.pPb_myOutput.root &
./lib/runGammaJets 2 0 true true user.jeouelle.2.4.30hi.calibcheck.signalonly.310.mc15_8TeV.423101.Slice2.pPb_myOutput.root &
./lib/runGammaJets 3 0 true true user.jeouelle.2.4.30hi.calibcheck.signalonly.310.mc15_8TeV.423102.Slice3.pPb_myOutput.root &
./lib/runGammaJets 4 0 true true user.jeouelle.2.4.30hi.calibcheck.signalonly.310.mc15_8TeV.423103.Slice4.pPb_myOutput.root &
./lib/runGammaJets 5 0 true true user.jeouelle.2.4.30hi.calibcheck.signalonly.310.mc15_8TeV.423104.Slice5.pPb_myOutput.root &
./lib/runGammaJets 6 0 true true user.jeouelle.2.4.30hi.calibcheck.signalonly.310.mc15_8TeV.423105.Slice6.pPb_myOutput.root &

wait

./lib/runGammaJetsHist  
