# gamma+jet pp sample with pPb data overlay
./lib/runEnergyScaleChecks 1 true user.jeouelle.2.4.30hi.calibcheck.200.mc15_8TeV.423100.Slice1.pPb_myOutput.root &
./lib/runEnergyScaleChecks 2 true user.jeouelle.2.4.30hi.calibcheck.200.mc15_8TeV.423101.Slice2.pPb_myOutput.root &
./lib/runEnergyScaleChecks 3 true user.jeouelle.2.4.30hi.calibcheck.200.mc15_8TeV.423102.Slice3.pPb_myOutput.root &
./lib/runEnergyScaleChecks 4 true user.jeouelle.2.4.30hi.calibcheck.200.mc15_8TeV.423103.Slice4.pPb_myOutput.root &
./lib/runEnergyScaleChecks 5 true user.jeouelle.2.4.30hi.calibcheck.200.mc15_8TeV.423104.Slice5.pPb_myOutput.root &
./lib/runEnergyScaleChecks 6 true user.jeouelle.2.4.30hi.calibcheck.200.mc15_8TeV.423105.Slice6.pPb_myOutput.root &

wait

./lib/runEnergyScaleChecks 1 false user.jeouelle.2.4.30hi.calibcheck.200.mc15_8TeV.423100.Slice1.Pbp_myOutput.root &
./lib/runEnergyScaleChecks 2 false user.jeouelle.2.4.30hi.calibcheck.200.mc15_8TeV.423101.Slice2.Pbp_myOutput.root &
./lib/runEnergyScaleChecks 3 false user.jeouelle.2.4.30hi.calibcheck.200.mc15_8TeV.423102.Slice3.Pbp_myOutput.root &
./lib/runEnergyScaleChecks 4 false user.jeouelle.2.4.30hi.calibcheck.200.mc15_8TeV.423103.Slice4.Pbp_myOutput.root &
./lib/runEnergyScaleChecks 5 false user.jeouelle.2.4.30hi.calibcheck.200.mc15_8TeV.423104.Slice5.Pbp_myOutput.root &
./lib/runEnergyScaleChecks 6 false user.jeouelle.2.4.30hi.calibcheck.200.mc15_8TeV.423105.Slice6.Pbp_myOutput.root &

wait

# Zmumu (data overlay) and Zee (signal only) samples for period A
./lib/runEnergyScaleChecks 0 true user.jeouelle.2.4.30hi.calibcheck.200.mc15_8TeV.361107.ZmumuJet.pPb_myOutput.root &
./lib/runEnergyScaleChecks 0 true user.jeouelle.2.4.30hi.calibcheck.200.mc15_8TeV.361106.ZeeJet.pPb_myOutput.root &

# ... and for period B
./lib/runEnergyScaleChecks 0 false user.jeouelle.2.4.30hi.calibcheck.200.mc15_8TeV.361107.ZmumuJet.Pbp_myOutput.root &
./lib/runEnergyScaleChecks 0 false user.jeouelle.2.4.30hi.calibcheck.200.mc15_8TeV.361106.ZeeJet.Pbp_myOutput.root &

wait

./lib/runEnergyScaleChecksHist  
