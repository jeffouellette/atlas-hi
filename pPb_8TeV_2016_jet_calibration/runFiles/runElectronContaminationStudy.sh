# Zee sample for period A
./lib/runElectronContaminationStudy true user.jeouelle.2.4.30hi.calibcheck.120.mc15_8TeV.361106.ZeeJet.pPb_myOutput.root &
# ... and for period B
./lib/runElectronContaminationStudy false user.jeouelle.2.4.30hi.calibcheck.120.mc15_8TeV.361106.ZeeJet.Pbp_myOutput.root &

wait

./lib/runElectronContaminationStudyHist
