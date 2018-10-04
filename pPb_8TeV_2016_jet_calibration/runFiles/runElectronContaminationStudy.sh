# Zee sample for period A
./bin/runElectronContaminationStudy true user.jeouelle.2.4.30hi.cabincheck.120.mc15_8TeV.361106.ZeeJet.pPb_myOutput.root &
# ... and for period B
./bin/runElectronContaminationStudy false user.jeouelle.2.4.30hi.cabincheck.120.mc15_8TeV.361106.ZeeJet.Pbp_myOutput.root &

wait

./bin/runElectronContaminationStudyHist
