version=400

# 2016 8TeV pPb main physics stream
./bin/runJetAnalysis data${version} 313063 0.03 false true &
./bin/runJetAnalysis data${version} 313067 1.24 false true &
./bin/runJetAnalysis data${version} 313100 9.66 false true &
./bin/runJetAnalysis data${version} 313107 11.92 false true & 
./bin/runJetAnalysis data${version} 313136 10.4 false true & 
./bin/runJetAnalysis data${version} 313187 3.67 false true &

wait

./bin/runJetAnalysis data${version} 313259 5.12 false true &
./bin/runJetAnalysis data${version} 313285 4.74 false true &
./bin/runJetAnalysis data${version} 313295 10.69 false true &
./bin/runJetAnalysis data${version} 313333 4.13 false true &
./bin/runJetAnalysis data${version} 313435 0.39 false true &
./bin/runJetAnalysis data${version} 313572 0.01 false false &

wait

./bin/runJetAnalysis data${version} 313574 1.33 false false &
./bin/runJetAnalysis data${version} 313575 7.54 false false &
./bin/runJetAnalysis data${version} 313603 8.69 false false & 
./bin/runJetAnalysis data${version} 313629 6.86 false false &
./bin/runJetAnalysis data${version} 313630 7.90 false false & 
./bin/runJetAnalysis data${version} 313688 7.96 false false &

wait

./bin/runJetAnalysis data${version} 313695 4.53 false false & 
./bin/runJetAnalysis data${version} 313833 5.11 false false &
./bin/runJetAnalysis data${version} 313878 2.16 false false & 
./bin/runJetAnalysis data${version} 313929 0.63 false false & 
./bin/runJetAnalysis data${version} 313935 10.96 false false &
./bin/runJetAnalysis data${version} 313984 2.40 false false &

wait

./bin/runJetAnalysis data${version} 314014 7.36 false false & 
./bin/runJetAnalysis data${version} 314077 10.19 false false &
./bin/runJetAnalysis data${version} 314105 6.50 false false &
./bin/runJetAnalysis data${version} 314112 10.49 false false &
./bin/runJetAnalysis data${version} 314157 9.83 false false & 
./bin/runJetAnalysis data${version} 314170 4.92 false false & 

wait

version=400

# gamma+jet pp sample with pPb data overlay
./bin/runJetAnalysis mc_${version} 1 0 true true user.jeouelle.2.4.30hi.calibcheck.${version}.mc15_8TeV.423100.Slice1.pPb_myOutput.root 9.5190E+03 0.000028838 306881 &
./bin/runJetAnalysis mc_${version} 2 0 true true user.jeouelle.2.4.30hi.calibcheck.${version}.mc15_8TeV.423101.Slice2.pPb_myOutput.root 7.2630E+02 0.000026548 311173 &
./bin/runJetAnalysis mc_${version} 3 0 true true user.jeouelle.2.4.30hi.calibcheck.${version}.mc15_8TeV.423102.Slice3.pPb_myOutput.root 1.8801E+02 0.000029816 313244 &
./bin/runJetAnalysis mc_${version} 4 0 true true user.jeouelle.2.4.30hi.calibcheck.${version}.mc15_8TeV.423103.Slice4.pPb_myOutput.root 5.0077E+01 0.000041955 311028 &
./bin/runJetAnalysis mc_${version} 5 0 true true user.jeouelle.2.4.30hi.calibcheck.${version}.mc15_8TeV.423104.Slice5.pPb_myOutput.root 2.8100E+00 0.000050105 312595 &
./bin/runJetAnalysis mc_${version} 6 0 true true user.jeouelle.2.4.30hi.calibcheck.${version}.mc15_8TeV.423105.Slice6.pPb_myOutput.root 1.2389E-01 0.000051776 298527 &

wait

./bin/runJetAnalysis mc_${version} 1 0 true false user.jeouelle.2.4.30hi.calibcheck.${version}.mc15_8TeV.423100.Slice1.Pbp_myOutput.root 9.5190E+03 0.00002885 644631 &
./bin/runJetAnalysis mc_${version} 2 0 true false user.jeouelle.2.4.30hi.calibcheck.${version}.mc15_8TeV.423101.Slice2.Pbp_myOutput.root 7.2630E+02 0.000026538 665323 &
./bin/runJetAnalysis mc_${version} 3 0 true false user.jeouelle.2.4.30hi.calibcheck.${version}.mc15_8TeV.423102.Slice3.Pbp_myOutput.root 1.8801E+02 0.00002988 643842 &
./bin/runJetAnalysis mc_${version} 4 0 true false user.jeouelle.2.4.30hi.calibcheck.${version}.mc15_8TeV.423103.Slice4.Pbp_myOutput.root 5.0077E+01 0.000041994 658656 &
./bin/runJetAnalysis mc_${version} 5 0 true false user.jeouelle.2.4.30hi.calibcheck.${version}.mc15_8TeV.423104.Slice5.Pbp_myOutput.root 2.8100E+00 0.000050105 665817 &
./bin/runJetAnalysis mc_${version} 6 0 true false user.jeouelle.2.4.30hi.calibcheck.${version}.mc15_8TeV.423105.Slice6.Pbp_myOutput.root 1.2388E-01 0.000051791 664131 &

wait

#version=310
#
## alternative validation signal-only gamma+jet sample (period A only!)
#./bin/runJetAnalysis mc_${version} 1 0 true true user.jeouelle.2.4.30hi.calibcheck.signalonly.${version}.mc15_8TeV.423100.Slice1.pPb_myOutput.root 9.5180E+03 0.00002901 49900 &
#./bin/runJetAnalysis mc_${version} 2 0 true true user.jeouelle.2.4.30hi.calibcheck.signalonly.${version}.mc15_8TeV.423101.Slice2.pPb_myOutput.root 7.26290E+02 0.000026722 49899 &
#./bin/runJetAnalysis mc_${version} 3 0 true true user.jeouelle.2.4.30hi.calibcheck.signalonly.${version}.mc15_8TeV.423102.Slice3.pPb_myOutput.root 1.88010E+02 0.000030018 50000 &
#./bin/runJetAnalysis mc_${version} 4 0 true true user.jeouelle.2.4.30hi.calibcheck.signalonly.${version}.mc15_8TeV.423103.Slice4.pPb_myOutput.root 5.0077E+01 0.000042496 49900 &
#./bin/runJetAnalysis mc_${version} 5 0 true true user.jeouelle.2.4.30hi.calibcheck.signalonly.${version}.mc15_8TeV.423104.Slice5.pPb_myOutput.root 2.810E+00 0.000051607 48000 &
#./bin/runJetAnalysis mc_${version} 6 0 true true user.jeouelle.2.4.30hi.calibcheck.signalonly.${version}.mc15_8TeV.423105.Slice6.pPb_myOutput.root 1.23880E-01 0.000052662 49900 &
#
#wait

rm ../rootFiles/JetAnalysis/outFile.root
hadd ../rootFiles/JetAnalysis/outFile.root ../rootFiles/JetAnalysis/dataSet_*.root

./bin/runJetAnalysisHist
