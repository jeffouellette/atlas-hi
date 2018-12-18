version=400

# 2016 8TeV pPb main physics stream
./bin/runZeeJets data${version} 313063 0.03 false true &
./bin/runZeeJets data${version} 313067 1.24 false true &
./bin/runZeeJets data${version} 313100 9.66 false true &
./bin/runZeeJets data${version} 313107 11.92 false true & 
./bin/runZeeJets data${version} 313136 10.4 false true & 
./bin/runZeeJets data${version} 313187 3.67 false true &

wait

./bin/runZeeJets data${version} 313259 5.12 false true &
./bin/runZeeJets data${version} 313285 4.74 false true &
./bin/runZeeJets data${version} 313295 10.69 false true &
./bin/runZeeJets data${version} 313333 4.13 false true &
./bin/runZeeJets data${version} 313435 0.39 false true &
./bin/runZeeJets data${version} 313572 0.01 false false &

wait

./bin/runZeeJets data${version} 313574 1.33 false false &
./bin/runZeeJets data${version} 313575 7.54 false false &
./bin/runZeeJets data${version} 313603 8.69 false false & 
./bin/runZeeJets data${version} 313629 6.86 false false &
./bin/runZeeJets data${version} 313630 7.90 false false & 
./bin/runZeeJets data${version} 313688 7.96 false false &

wait

./bin/runZeeJets data${version} 313695 4.53 false false & 
./bin/runZeeJets data${version} 313833 5.11 false false &
./bin/runZeeJets data${version} 313878 2.16 false false & 
./bin/runZeeJets data${version} 313929 0.63 false false & 
./bin/runZeeJets data${version} 313935 10.96 false false &
./bin/runZeeJets data${version} 313984 2.40 false false &

wait

./bin/runZeeJets data${version} 314014 7.36 false false & 
./bin/runZeeJets data${version} 314077 10.19 false false &
./bin/runZeeJets data${version} 314105 6.50 false false &
./bin/runZeeJets data${version} 314112 10.49 false false &
./bin/runZeeJets data${version} 314157 9.83 false false & 
./bin/runZeeJets data${version} 314170 4.92 false false & 

wait

version=400

# Zee (signal only) sample for period A
./bin/runZeeJets mc_${version} 0 0 true true user.jeouelle.2.4.30hi.calibcheck.${version}.mc15_8TeV.361106.ZeeJet.pPb_myOutput.root 1.1360E-03 1 167761 &

# ... and for period B
./bin/runZeeJets mc_${version} 0 0 true false user.jeouelle.2.4.30hi.calibcheck.${version}.mc15_8TeV.361106.ZeeJet.Pbp_myOutput.root 1.1360E-03 1 302313 &

wait

#version=420
#
## Alternative signal only Zee sample for period A
#./bin/runZeeJets mc_${version} 1 0 true true user.jeouelle.2.4.30hi.cabincheck.signalonly.420.mc15_8TeV.361106.ZeeJet1.pPb_myOutput.root &
#./bin/runZeeJets mc_${version} 2 0 true true user.jeouelle.2.4.30hi.cabincheck.signalonly.420.mc15_8TeV.361106.ZeeJet2.pPb_myOutput.root &
#./bin/runZeeJets mc_${version} 3 0 true true user.jeouelle.2.4.30hi.cabincheck.signalonly.420.mc15_8TeV.361106.ZeeJet3.pPb_myOutput.root &
#./bin/runZeeJets mc_${version} 4 0 true true user.jeouelle.2.4.30hi.cabincheck.signalonly.420.mc15_8TeV.361106.ZeeJet4.pPb_myOutput.root &
#./bin/runZeeJets mc_${version} 5 0 true true user.jeouelle.2.4.30hi.cabincheck.signalonly.420.mc15_8TeV.361106.ZeeJet5.pPb_myOutput.root &
#
#wait
#
## ... and for period B
#./bin/runZeeJets mc_${version} 1 0 true false user.jeouelle.2.4.30hi.cabincheck.signalonly.420.mc15_8TeV.361106.ZeeJet1.Pbp_myOutput.root &
#./bin/runZeeJets mc_${version} 2 0 true false user.jeouelle.2.4.30hi.cabincheck.signalonly.420.mc15_8TeV.361106.ZeeJet2.Pbp_myOutput.root &
#./bin/runZeeJets mc_${version} 3 0 true false user.jeouelle.2.4.30hi.cabincheck.signalonly.420.mc15_8TeV.361106.ZeeJet3.Pbp_myOutput.root &
#./bin/runZeeJets mc_${version} 4 0 true false user.jeouelle.2.4.30hi.cabincheck.signalonly.420.mc15_8TeV.361106.ZeeJet4.Pbp_myOutput.root &
#./bin/runZeeJets mc_${version} 5 0 true false user.jeouelle.2.4.30hi.cabincheck.signalonly.420.mc15_8TeV.361106.ZeeJet5.Pbp_myOutput.root &
#
#wait

rm ../rootFiles/ZeeJets/outFile.root
hadd ../rootFiles/ZeeJets/outFile.root ../rootFiles/ZeeJets/dataSet_*.root

./bin/runZeeJetsHist
