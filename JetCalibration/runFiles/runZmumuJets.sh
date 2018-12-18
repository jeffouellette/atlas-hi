version=400

# 2016 8TeV pPb main physics stream
./bin/runZmumuJets data${version} 313063 0.03 false true &
./bin/runZmumuJets data${version} 313067 1.24 false true &
./bin/runZmumuJets data${version} 313100 9.66 false true &
./bin/runZmumuJets data${version} 313107 11.92 false true & 
./bin/runZmumuJets data${version} 313136 10.4 false true & 
./bin/runZmumuJets data${version} 313187 3.67 false true &

wait

./bin/runZmumuJets data${version} 313259 5.12 false true &
./bin/runZmumuJets data${version} 313285 4.74 false true &
./bin/runZmumuJets data${version} 313295 10.69 false true &
./bin/runZmumuJets data${version} 313333 4.13 false true &
./bin/runZmumuJets data${version} 313435 0.39 false true &
./bin/runZmumuJets data${version} 313572 0.01 false false &

wait

./bin/runZmumuJets data${version} 313574 1.33 false false &
./bin/runZmumuJets data${version} 313575 7.54 false false &
./bin/runZmumuJets data${version} 313603 8.69 false false & 
./bin/runZmumuJets data${version} 313629 6.86 false false &
./bin/runZmumuJets data${version} 313630 7.90 false false & 
./bin/runZmumuJets data${version} 313688 7.96 false false &

wait

./bin/runZmumuJets data${version} 313695 4.53 false false & 
./bin/runZmumuJets data${version} 313833 5.11 false false &
./bin/runZmumuJets data${version} 313878 2.16 false false & 
./bin/runZmumuJets data${version} 313929 0.63 false false & 
./bin/runZmumuJets data${version} 313935 10.96 false false &
./bin/runZmumuJets data${version} 313984 2.40 false false &

wait

./bin/runZmumuJets data${version} 314014 7.36 false false & 
./bin/runZmumuJets data${version} 314077 10.19 false false &
./bin/runZmumuJets data${version} 314105 6.50 false false &
./bin/runZmumuJets data${version} 314112 10.49 false false &
./bin/runZmumuJets data${version} 314157 9.83 false false & 
./bin/runZmumuJets data${version} 314170 4.92 false false & 

wait

version=400

# Zmumu (data overlay) samples for period A
./bin/runZmumuJets mc_${version} 0 0 true true user.jeouelle.2.4.30hi.calibcheck.${version}.mc15_8TeV.361107.ZmumuJet.pPb_myOutput.root 1.1360E-03 1 131944 &

# ... and for period B
./bin/runZmumuJets mc_${version} 0 0 true false user.jeouelle.2.4.30hi.calibcheck.${version}.mc15_8TeV.361107.ZmumuJet.Pbp_myOutput.root 1.1360E-03 1 243320 &

wait

rm ../rootFiles/ZmumuJets/outFile.root
hadd ../rootFiles/ZmumuJets/outFile.root ../rootFiles/ZmumuJets/dataSet_*.root

./bin/runZmumuJetsHist
