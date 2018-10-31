#version=230
#
#./bin/runEnergyScaleChecks mc_${version} 2 true user.jeouelle.2.4.30hi.calibcheck.signalonly.${version}.mc15_8TeV.420012.jetjet.JZ2R04.pPb_myOutput.root 1.291E3 0.0056462 3997692 &
#./bin/runEnergyScaleChecks jetjet_valid_${version} 2 true user.jeouelle.2.4.30hi.calibcheck.signalonly.valid.${version}.mc15_8TeV.420012.jetjet.JZ2R04.pPb_myOutput.root 1.292E3 0.0056933 9500 &
#./bin/runEnergyScaleChecks jetjet_valid_${version} 2 false user.jeouelle.2.4.30hi.calibcheck.signalonly.valid.${version}.mc15_8TeV.420012.jetjet.JZ2R04.Pbp_myOutput.root 1.292E3 0.0056911 10000 &

version=200

# gamma+jet pp sample with pPb data overlay
./bin/runEnergyScaleChecks mc_${version} 1 true user.jeouelle.2.4.30hi.calibcheck.${version}.mc15_8TeV.423100.Slice1.pPb_myOutput.root &
./bin/runEnergyScaleChecks mc_${version} 2 true user.jeouelle.2.4.30hi.calibcheck.${version}.mc15_8TeV.423101.Slice2.pPb_myOutput.root &
./bin/runEnergyScaleChecks mc_${version} 3 true user.jeouelle.2.4.30hi.calibcheck.${version}.mc15_8TeV.423102.Slice3.pPb_myOutput.root &
./bin/runEnergyScaleChecks mc_${version} 4 true user.jeouelle.2.4.30hi.calibcheck.${version}.mc15_8TeV.423103.Slice4.pPb_myOutput.root &
./bin/runEnergyScaleChecks mc_${version} 5 true user.jeouelle.2.4.30hi.calibcheck.${version}.mc15_8TeV.423104.Slice5.pPb_myOutput.root &
./bin/runEnergyScaleChecks mc_${version} 6 true user.jeouelle.2.4.30hi.calibcheck.${version}.mc15_8TeV.423105.Slice6.pPb_myOutput.root &

wait

./bin/runEnergyScaleChecks mc_${version} 1 false user.jeouelle.2.4.30hi.calibcheck.${version}.mc15_8TeV.423100.Slice1.Pbp_myOutput.root &
./bin/runEnergyScaleChecks mc_${version} 2 false user.jeouelle.2.4.30hi.calibcheck.${version}.mc15_8TeV.423101.Slice2.Pbp_myOutput.root &
./bin/runEnergyScaleChecks mc_${version} 3 false user.jeouelle.2.4.30hi.calibcheck.${version}.mc15_8TeV.423102.Slice3.Pbp_myOutput.root &
./bin/runEnergyScaleChecks mc_${version} 4 false user.jeouelle.2.4.30hi.calibcheck.${version}.mc15_8TeV.423103.Slice4.Pbp_myOutput.root &
./bin/runEnergyScaleChecks mc_${version} 5 false user.jeouelle.2.4.30hi.calibcheck.${version}.mc15_8TeV.423104.Slice5.Pbp_myOutput.root &
./bin/runEnergyScaleChecks mc_${version} 6 false user.jeouelle.2.4.30hi.calibcheck.${version}.mc15_8TeV.423105.Slice6.Pbp_myOutput.root &

wait

version=310

# alternative validation signal-only gamma+jet sample (period A only!)
./bin/runEnergyScaleChecks mc_${version} 1 true user.jeouelle.2.4.30hi.calibcheck.signalonly.${version}.mc15_8TeV.423100.Slice1.pPb_myOutput.root &
./bin/runEnergyScaleChecks mc_${version} 2 true user.jeouelle.2.4.30hi.calibcheck.signalonly.${version}.mc15_8TeV.423101.Slice2.pPb_myOutput.root &
./bin/runEnergyScaleChecks mc_${version} 3 true user.jeouelle.2.4.30hi.calibcheck.signalonly.${version}.mc15_8TeV.423102.Slice3.pPb_myOutput.root &
./bin/runEnergyScaleChecks mc_${version} 4 true user.jeouelle.2.4.30hi.calibcheck.signalonly.${version}.mc15_8TeV.423103.Slice4.pPb_myOutput.root &
./bin/runEnergyScaleChecks mc_${version} 5 true user.jeouelle.2.4.30hi.calibcheck.signalonly.${version}.mc15_8TeV.423104.Slice5.pPb_myOutput.root &
./bin/runEnergyScaleChecks mc_${version} 6 true user.jeouelle.2.4.30hi.calibcheck.signalonly.${version}.mc15_8TeV.423105.Slice6.pPb_myOutput.root &

wait

version=200

# Zmumu (data overlay) and Zee (signal only) samples for period A
./bin/runEnergyScaleChecks mc_${version} 0 true user.jeouelle.2.4.30hi.calibcheck.${version}.mc15_8TeV.361107.ZmumuJet.pPb_myOutput.root &
./bin/runEnergyScaleChecks mc_${version} 0 true user.jeouelle.2.4.30hi.calibcheck.${version}.mc15_8TeV.361106.ZeeJet.pPb_myOutput.root &

# ... and for period B
./bin/runEnergyScaleChecks mc_${version} 0 false user.jeouelle.2.4.30hi.calibcheck.${version}.mc15_8TeV.361107.ZmumuJet.Pbp_myOutput.root &
./bin/runEnergyScaleChecks mc_${version} 0 false user.jeouelle.2.4.30hi.calibcheck.${version}.mc15_8TeV.361106.ZeeJet.Pbp_myOutput.root &

wait

./bin/runEnergyScaleChecksHist
