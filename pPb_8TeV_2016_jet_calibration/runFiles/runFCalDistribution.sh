# 2016 8TeV pPb main physics stream
./bin/runFCalDistribution 313063 0.03 false true &
./bin/runFCalDistribution 313067 1.24 false true &
./bin/runFCalDistribution 313100 9.66 false true &
./bin/runFCalDistribution 313107 11.92 false true & 
./bin/runFCalDistribution 313136 10.4 false true & 
./bin/runFCalDistribution 313187 3.67 false true &

wait

./bin/runFCalDistribution 313259 5.12 false true &
./bin/runFCalDistribution 313285 4.74 false true &
./bin/runFCalDistribution 313295 10.69 false true &
./bin/runFCalDistribution 313333 4.13 false true &
./bin/runFCalDistribution 313435 0.39 false true &
./bin/runFCalDistribution 313572 0.01 false false &

wait

./bin/runFCalDistribution 313574 1.33 false false &
./bin/runFCalDistribution 313575 7.54 false false &
./bin/runFCalDistribution 313603 8.69 false false & 
./bin/runFCalDistribution 313629 6.86 false false &
./bin/runFCalDistribution 313630 7.90 false false & 
./bin/runFCalDistribution 313688 7.96 false false &

wait

./bin/runFCalDistribution 313695 4.53 false false & 
./bin/runFCalDistribution 313833 5.11 false false &
./bin/runFCalDistribution 313878 2.16 false false & 
./bin/runFCalDistribution 313929 0.63 false false & 
./bin/runFCalDistribution 313935 10.96 false false &
./bin/runFCalDistribution 313984 2.40 false false &

wait

./bin/runFCalDistribution 314014 7.36 false false & 
./bin/runFCalDistribution 314077 10.19 false false &
./bin/runFCalDistribution 314105 6.50 false false &
./bin/runFCalDistribution 314112 10.49 false false &
./bin/runFCalDistribution 314157 9.83 false false & 
./bin/runFCalDistribution 314170 4.92 false false & 

wait

# gamma+jet pp sample with pPb data overlay
./bin/runFCalDistribution 1 0 true true user.jeouelle.2.4.30hi.cabincheck.200.mc15_8TeV.423100.Slice1.pPb_myOutput.root &
./bin/runFCalDistribution 2 0 true true user.jeouelle.2.4.30hi.cabincheck.200.mc15_8TeV.423101.Slice2.pPb_myOutput.root &
./bin/runFCalDistribution 3 0 true true user.jeouelle.2.4.30hi.cabincheck.200.mc15_8TeV.423102.Slice3.pPb_myOutput.root &
./bin/runFCalDistribution 4 0 true true user.jeouelle.2.4.30hi.cabincheck.200.mc15_8TeV.423103.Slice4.pPb_myOutput.root &
./bin/runFCalDistribution 5 0 true true user.jeouelle.2.4.30hi.cabincheck.200.mc15_8TeV.423104.Slice5.pPb_myOutput.root &
./bin/runFCalDistribution 6 0 true true user.jeouelle.2.4.30hi.cabincheck.200.mc15_8TeV.423105.Slice6.pPb_myOutput.root &

wait

./bin/runFCalDistribution 1 0 true false user.jeouelle.2.4.30hi.cabincheck.200.mc15_8TeV.423100.Slice1.Pbp_myOutput.root &
./bin/runFCalDistribution 2 0 true false user.jeouelle.2.4.30hi.cabincheck.200.mc15_8TeV.423101.Slice2.Pbp_myOutput.root &
./bin/runFCalDistribution 3 0 true false user.jeouelle.2.4.30hi.cabincheck.200.mc15_8TeV.423102.Slice3.Pbp_myOutput.root &
./bin/runFCalDistribution 4 0 true false user.jeouelle.2.4.30hi.cabincheck.200.mc15_8TeV.423103.Slice4.Pbp_myOutput.root &
./bin/runFCalDistribution 5 0 true false user.jeouelle.2.4.30hi.cabincheck.200.mc15_8TeV.423104.Slice5.Pbp_myOutput.root &
./bin/runFCalDistribution 6 0 true false user.jeouelle.2.4.30hi.cabincheck.200.mc15_8TeV.423105.Slice6.Pbp_myOutput.root &

wait

# Zmumu (data overlay) and Zee (signal only) samples for period A
./bin/runFCalDistribution 0 0 true true user.jeouelle.2.4.30hi.cabincheck.200.mc15_8TeV.361107.ZmumuJet.pPb_myOutput.root &
./bin/runFCalDistribution 0 0 true true user.jeouelle.2.4.30hi.cabincheck.200.mc15_8TeV.361106.ZeeJet.pPb_myOutput.root &

# ... and for period B
./bin/runFCalDistribution 0 0 true false user.jeouelle.2.4.30hi.cabincheck.200.mc15_8TeV.361107.ZmumuJet.Pbp_myOutput.root &
./bin/runFCalDistribution 0 0 true false user.jeouelle.2.4.30hi.cabincheck.200.mc15_8TeV.361106.ZeeJet.Pbp_myOutput.root &

wait

./bin/runFCalDistributionHist
