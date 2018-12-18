version=400

# 2016 8TeV pPb main physics stream
./bin/runFCalDistribution data${version} 313063 0.03 false true &
./bin/runFCalDistribution data${version} 313067 1.24 false true &
./bin/runFCalDistribution data${version} 313100 9.66 false true &
./bin/runFCalDistribution data${version} 313107 11.92 false true & 
./bin/runFCalDistribution data${version} 313136 10.4 false true & 
./bin/runFCalDistribution data${version} 313187 3.67 false true &

wait

./bin/runFCalDistribution data${version} 313259 5.12 false true &
./bin/runFCalDistribution data${version} 313285 4.74 false true &
./bin/runFCalDistribution data${version} 313295 10.69 false true &
./bin/runFCalDistribution data${version} 313333 4.13 false true &
./bin/runFCalDistribution data${version} 313435 0.39 false true &
./bin/runFCalDistribution data${version} 313572 0.01 false false &

wait

./bin/runFCalDistribution data${version} 313574 1.33 false false &
./bin/runFCalDistribution data${version} 313575 7.54 false false &
./bin/runFCalDistribution data${version} 313603 8.69 false false & 
./bin/runFCalDistribution data${version} 313629 6.86 false false &
./bin/runFCalDistribution data${version} 313630 7.90 false false & 
./bin/runFCalDistribution data${version} 313688 7.96 false false &

wait

./bin/runFCalDistribution data${version} 313695 4.53 false false & 
./bin/runFCalDistribution data${version} 313833 5.11 false false &
./bin/runFCalDistribution data${version} 313878 2.16 false false & 
./bin/runFCalDistribution data${version} 313929 0.63 false false & 
./bin/runFCalDistribution data${version} 313935 10.96 false false &
./bin/runFCalDistribution data${version} 313984 2.40 false false &

wait

./bin/runFCalDistribution data${version} 314014 7.36 false false & 
./bin/runFCalDistribution data${version} 314077 10.19 false false &
./bin/runFCalDistribution data${version} 314105 6.50 false false &
./bin/runFCalDistribution data${version} 314112 10.49 false false &
./bin/runFCalDistribution data${version} 314157 9.83 false false & 
./bin/runFCalDistribution data${version} 314170 4.92 false false & 

wait

version=400

# gamma+jet pp sample with pPb data overlay
./bin/runFCalDistribution mc_${version} 1 0 true true user.jeouelle.2.4.30hi.calibcheck.${version}.mc15_8TeV.423100.Slice1.pPb_myOutput.root 9.5190E+03 0.000028838 306881 &
./bin/runFCalDistribution mc_${version} 2 0 true true user.jeouelle.2.4.30hi.calibcheck.${version}.mc15_8TeV.423101.Slice2.pPb_myOutput.root 7.2630E+02 0.000026548 311173 &
./bin/runFCalDistribution mc_${version} 3 0 true true user.jeouelle.2.4.30hi.calibcheck.${version}.mc15_8TeV.423102.Slice3.pPb_myOutput.root 1.8801E+02 0.000029816 313244 &
./bin/runFCalDistribution mc_${version} 4 0 true true user.jeouelle.2.4.30hi.calibcheck.${version}.mc15_8TeV.423103.Slice4.pPb_myOutput.root 5.0077E+01 0.000041955 311028 &
./bin/runFCalDistribution mc_${version} 5 0 true true user.jeouelle.2.4.30hi.calibcheck.${version}.mc15_8TeV.423104.Slice5.pPb_myOutput.root 2.8100E+00 0.000050105 312595 &
./bin/runFCalDistribution mc_${version} 6 0 true true user.jeouelle.2.4.30hi.calibcheck.${version}.mc15_8TeV.423105.Slice6.pPb_myOutput.root 1.2389E-01 0.000051776 298527 &

wait

./bin/runFCalDistribution mc_${version} 1 0 true false user.jeouelle.2.4.30hi.calibcheck.${version}.mc15_8TeV.423100.Slice1.Pbp_myOutput.root 9.5190E+03 0.00002885 644631 &
./bin/runFCalDistribution mc_${version} 2 0 true false user.jeouelle.2.4.30hi.calibcheck.${version}.mc15_8TeV.423101.Slice2.Pbp_myOutput.root 7.2630E+02 0.000026538 665323 &
./bin/runFCalDistribution mc_${version} 3 0 true false user.jeouelle.2.4.30hi.calibcheck.${version}.mc15_8TeV.423102.Slice3.Pbp_myOutput.root 1.8801E+02 0.00002988 643842 &
./bin/runFCalDistribution mc_${version} 4 0 true false user.jeouelle.2.4.30hi.calibcheck.${version}.mc15_8TeV.423103.Slice4.Pbp_myOutput.root 5.0077E+01 0.000041994 658656 &
./bin/runFCalDistribution mc_${version} 5 0 true false user.jeouelle.2.4.30hi.calibcheck.${version}.mc15_8TeV.423104.Slice5.Pbp_myOutput.root 2.8100E+00 0.000050105 665817 &
./bin/runFCalDistribution mc_${version} 6 0 true false user.jeouelle.2.4.30hi.calibcheck.${version}.mc15_8TeV.423105.Slice6.Pbp_myOutput.root 1.2388E-01 0.000051791 664131 &

wait

## Zmumu (data overlay) samples for period A
#./bin/runFCalDistribution mc_${version} 0 0 true true user.jeouelle.2.4.30hi.calibcheck.${version}.mc15_8TeV.361107.ZmumuJet.pPb_myOutput.root 1.1360E-03 1 131944 &
## ... and for period B
#./bin/runFCalDistribution mc_${version} 0 0 true false user.jeouelle.2.4.30hi.calibcheck.${version}.mc15_8TeV.361107.ZmumuJet.Pbp_myOutput.root 1.1360E-03 1 243320 &
#
## Zee (signal only) sample for period A
#./bin/runFCalDistribution mc_${version} 0 true true user.jeouelle.2.4.30hi.calibcheck.${version}.mc15_8TeV.361106.ZeeJet.pPb_myOutput.root 1.1360E-03 1 167761 &
## ... and for period B
#./bin/runFCalDistribution mc_${version} 0 true false user.jeouelle.2.4.30hi.calibcheck.${version}.mc15_8TeV.361106.ZeeJet.Pbp_myOutput.root 1.1360E-03 1 302313 &
#
#wait

rm ../rootFiles/FCalDistribution/outFile.root
hadd ../rootFiles/FCalDistribution/outFile.root ../rootFiles/FCalDistribution/dataSet_*.root

./bin/runFCalDistributionHist
