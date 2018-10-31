#version=230
#
#./bin/runIntercalibration mc_${version} 2 true 60 160 user.jeouelle.2.4.30hi.calibcheck.signalonly.${version}.mc15_8TeV.420012.jetjet.JZ2R04.pPb_myOutput.root false 1.291E3 0.0056462 3997692 &

version=410
./bin/runIntercalibration jetjet_valid_${version} 2 true 60 160 user.jeouelle.2.4.30hi.calibcheck.valid.${version}.mc15_8TeV.420012.jetjet.JZ2R04.pPb_myOutput.root false 1 1 9875 &
./bin/runIntercalibration jetjet_valid_${version} 2 false 60 160 user.jeouelle.2.4.30hi.calibcheck.valid.${version}.mc15_8TeV.420012.jetjet.JZ2R04.Pbp_myOutput.root false 1 1 9835 &

wait

#version=270
#
## gamma+jet pp sample with pPb data overlay
#./bin/runIntercalibration mc_${version} 1 true 17 35 user.jeouelle.2.4.30hi.calibcheck.${version}.mc15_8TeV.423100.Slice1.pPb_myOutput.root true 9.5190E+03 0.000028838 306881 &
#./bin/runIntercalibration mc_${version} 2 true 35 50 user.jeouelle.2.4.30hi.calibcheck.${version}.mc15_8TeV.423101.Slice2.pPb_myOutput.root true 7.2630E+02 0.000026548 311173 &
#./bin/runIntercalibration mc_${version} 3 true 50 70 user.jeouelle.2.4.30hi.calibcheck.${version}.mc15_8TeV.423102.Slice3.pPb_myOutput.root true 1.8801E+02 0.000029816 313244 &
#./bin/runIntercalibration mc_${version} 4 true 70 140 user.jeouelle.2.4.30hi.calibcheck.${version}.mc15_8TeV.423103.Slice4.pPb_myOutput.root true 5.0077E+01 0.000041955 311028 &
#./bin/runIntercalibration mc_${version} 5 true 140 280 user.jeouelle.2.4.30hi.calibcheck.${version}.mc15_8TeV.423104.Slice5.pPb_myOutput.root true 2.8100E+00 0.000050105 312595 &
#./bin/runIntercalibration mc_${version} 6 true 280 500 user.jeouelle.2.4.30hi.calibcheck.${version}.mc15_8TeV.423105.Slice6.pPb_myOutput.root true 1.2389E-01 0.000051776 298527 &
#
#wait
#
#./bin/runIntercalibration mc_${version} 1 false 17 35 user.jeouelle.2.4.30hi.calibcheck.${version}.mc15_8TeV.423100.Slice1.Pbp_myOutput.root true 9.5190E+03 0.00002885 644631 &
#./bin/runIntercalibration mc_${version} 2 false 35 50 user.jeouelle.2.4.30hi.calibcheck.${version}.mc15_8TeV.423101.Slice2.Pbp_myOutput.root true 7.2630E+02 0.000026538 665323 &
#./bin/runIntercalibration mc_${version} 3 false 50 70 user.jeouelle.2.4.30hi.calibcheck.${version}.mc15_8TeV.423102.Slice3.Pbp_myOutput.root true 1.8801E+02 0.00002988 643842 &
#./bin/runIntercalibration mc_${version} 4 false 70 140 user.jeouelle.2.4.30hi.calibcheck.${version}.mc15_8TeV.423103.Slice4.Pbp_myOutput.root true 5.0077E+01 0.000041994 658656 &
#./bin/runIntercalibration mc_${version} 5 false 140 280 user.jeouelle.2.4.30hi.calibcheck.${version}.mc15_8TeV.423104.Slice5.Pbp_myOutput.root true 2.8100E+00 0.000050105 665817 &
#./bin/runIntercalibration mc_${version} 6 false 280 500 user.jeouelle.2.4.30hi.calibcheck.${version}.mc15_8TeV.423105.Slice6.Pbp_myOutput.root true 1.2388E-01 0.000051791 664131 &
#
#wait
#
#version=310
#
## alternative validation signal-only gamma+jet sample (period A only!)
#./bin/runIntercalibration mc_${version} 1 true 17 35 user.jeouelle.2.4.30hi.calibcheck.signalonly.${version}.mc15_8TeV.423100.Slice1.pPb_myOutput.root true 9.5180E+03 0.00002901 49900 &
#./bin/runIntercalibration mc_${version} 2 true 35 50 user.jeouelle.2.4.30hi.calibcheck.signalonly.${version}.mc15_8TeV.423101.Slice2.pPb_myOutput.root true 7.26290E+02 0.000026722 49899 &
#./bin/runIntercalibration mc_${version} 3 true 50 70 user.jeouelle.2.4.30hi.calibcheck.signalonly.${version}.mc15_8TeV.423102.Slice3.pPb_myOutput.root true 1.88010E+02 0.000030018 50000 &
#./bin/runIntercalibration mc_${version} 4 true 70 140 user.jeouelle.2.4.30hi.calibcheck.signalonly.${version}.mc15_8TeV.423103.Slice4.pPb_myOutput.root true 5.0077E+01 0.000042496 49900 &
#./bin/runIntercalibration mc_${version} 5 true 140 280 user.jeouelle.2.4.30hi.calibcheck.signalonly.${version}.mc15_8TeV.423104.Slice5.pPb_myOutput.root true 2.810E+00 0.000051607 48000 &
#./bin/runIntercalibration mc_${version} 6 true 280 500 user.jeouelle.2.4.30hi.calibcheck.signalonly.${version}.mc15_8TeV.423105.Slice6.pPb_myOutput.root true 1.23880E-01 0.000052662 49900 &
#
#wait

./bin/runIntercalibrationHist
