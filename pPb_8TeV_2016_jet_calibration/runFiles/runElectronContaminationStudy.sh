# Zee sample for period A
root -l -b -q '../AnalysisCode/electronContaminationStudy.C (-6, true, "user.jeouelle.2.4.30hi.calibcheck.900.mc15_8TeV.361106.ZeeJet.pPb_myOutput.root")'&
# ... and for period B
root -l -b -q '../AnalysisCode/electronContaminationStudy.C (-6, false, "user.jeouelle.2.4.30hi.calibcheck.900.mc15_8TeV.361106.ZeeJet.Pbp_myOutput.root")'&

wait

root -l -b -q '../AnalysisCode/electronContaminationStudyHist.C ()'
