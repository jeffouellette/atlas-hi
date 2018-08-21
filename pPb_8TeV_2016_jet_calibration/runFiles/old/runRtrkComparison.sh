# 2016 8TeV pPb main physics stream
root -l -b -q '../AnalysisCode/RtrkComparison.C (313063, 0.03, false, true)'&
root -l -b -q '../AnalysisCode/RtrkComparison.C (313067, 1.24, false, true)'&
root -l -b -q '../AnalysisCode/RtrkComparison.C (313100, 9.66, false, true)'&
root -l -b -q '../AnalysisCode/RtrkComparison.C (313107, 11.92, false, true)'& 
root -l -b -q '../AnalysisCode/RtrkComparison.C (313136, 10.4, false, true)'& 
root -l -b -q '../AnalysisCode/RtrkComparison.C (313187, 3.67, false, true)'&

wait

root -l -b -q '../AnalysisCode/RtrkComparison.C (313259, 5.12, false, true)'&
root -l -b -q '../AnalysisCode/RtrkComparison.C (313285, 4.74, false, true)'&
root -l -b -q '../AnalysisCode/RtrkComparison.C (313295, 10.69, false, true)'&
root -l -b -q '../AnalysisCode/RtrkComparison.C (313333, 4.13, false, true)'&
root -l -b -q '../AnalysisCode/RtrkComparison.C (313435, 0.39, false, true)'&
root -l -b -q '../AnalysisCode/RtrkComparison.C (313572, 0.01, false, false)'&

wait

root -l -b -q '../AnalysisCode/RtrkComparison.C (313574, 1.33, false, false)'&
root -l -b -q '../AnalysisCode/RtrkComparison.C (313575, 7.54, false, false)'&
root -l -b -q '../AnalysisCode/RtrkComparison.C (313603, 8.69, false, false)'& 
root -l -b -q '../AnalysisCode/RtrkComparison.C (313629, 6.86, false, false)'&
root -l -b -q '../AnalysisCode/RtrkComparison.C (313630, 7.90, false, false)'& 
root -l -b -q '../AnalysisCode/RtrkComparison.C (313688, 7.96, false, false)'&

wait

root -l -b -q '../AnalysisCode/RtrkComparison.C (313695, 4.53, false, false)'& 
root -l -b -q '../AnalysisCode/RtrkComparison.C (313833, 5.11, false, false)'&
root -l -b -q '../AnalysisCode/RtrkComparison.C (313878, 2.16, false, false)'& 
root -l -b -q '../AnalysisCode/RtrkComparison.C (313929, 0.63, false, false)'& 
root -l -b -q '../AnalysisCode/RtrkComparison.C (313935, 10.96, false, false)'&
root -l -b -q '../AnalysisCode/RtrkComparison.C (313984, 2.40, false, false)'&

wait

root -l -b -q '../AnalysisCode/RtrkComparison.C (314014, 7.36, false, false)'& 
root -l -b -q '../AnalysisCode/RtrkComparison.C (314077, 10.19, false, false)'&
root -l -b -q '../AnalysisCode/RtrkComparison.C (314105, 6.50, false, false)'&
root -l -b -q '../AnalysisCode/RtrkComparison.C (314112, 10.49, false, false)'&
root -l -b -q '../AnalysisCode/RtrkComparison.C (314157, 9.83, false, false)'& 
root -l -b -q '../AnalysisCode/RtrkComparison.C (314170, 4.92, false, false)'& 

wait

# gamma+jet pp sample with pPb data overlay
root -l -b -q '../AnalysisCode/RtrkComparison.C (1, 0, true, true, "user.jeouelle.2.4.30hi.calibcheck.200.mc15_8TeV.423100.Slice1.pPb_myOutput.root")'&
root -l -b -q '../AnalysisCode/RtrkComparison.C (2, 0, true, true, "user.jeouelle.2.4.30hi.calibcheck.200.mc15_8TeV.423101.Slice2.pPb_myOutput.root")'&
root -l -b -q '../AnalysisCode/RtrkComparison.C (3, 0, true, true, "user.jeouelle.2.4.30hi.calibcheck.200.mc15_8TeV.423102.Slice3.pPb_myOutput.root")'&
root -l -b -q '../AnalysisCode/RtrkComparison.C (4, 0, true, true, "user.jeouelle.2.4.30hi.calibcheck.200.mc15_8TeV.423103.Slice4.pPb_myOutput.root")'&
root -l -b -q '../AnalysisCode/RtrkComparison.C (5, 0, true, true, "user.jeouelle.2.4.30hi.calibcheck.200.mc15_8TeV.423104.Slice5.pPb_myOutput.root")'&
root -l -b -q '../AnalysisCode/RtrkComparison.C (6, 0, true, true, "user.jeouelle.2.4.30hi.calibcheck.200.mc15_8TeV.423105.Slice6.pPb_myOutput.root")'&

wait

root -l -b -q '../AnalysisCode/RtrkComparison.C (1, 0, true, false, "user.jeouelle.2.4.30hi.calibcheck.200.mc15_8TeV.423100.Slice1.Pbp_myOutput.root")'&
root -l -b -q '../AnalysisCode/RtrkComparison.C (2, 0, true, false, "user.jeouelle.2.4.30hi.calibcheck.200.mc15_8TeV.423101.Slice2.Pbp_myOutput.root")'&
root -l -b -q '../AnalysisCode/RtrkComparison.C (3, 0, true, false, "user.jeouelle.2.4.30hi.calibcheck.200.mc15_8TeV.423102.Slice3.Pbp_myOutput.root")'&
root -l -b -q '../AnalysisCode/RtrkComparison.C (4, 0, true, false, "user.jeouelle.2.4.30hi.calibcheck.200.mc15_8TeV.423103.Slice4.Pbp_myOutput.root")'&
root -l -b -q '../AnalysisCode/RtrkComparison.C (5, 0, true, false, "user.jeouelle.2.4.30hi.calibcheck.200.mc15_8TeV.423104.Slice5.Pbp_myOutput.root")'&
root -l -b -q '../AnalysisCode/RtrkComparison.C (6, 0, true, false, "user.jeouelle.2.4.30hi.calibcheck.200.mc15_8TeV.423105.Slice6.Pbp_myOutput.root")'&

wait

# Zmumu (data overlay) and Zee (signal only) samples for period A
root -l -b -q '../AnalysisCode/RtrkComparison.C (0, 0, true, true, "user.jeouelle.2.4.30hi.calibcheck.200.mc15_8TeV.361107.ZmumuJet.pPb_myOutput.root")'&
root -l -b -q '../AnalysisCode/RtrkComparison.C (-6, 0, true, true, "user.jeouelle.2.4.30hi.calibcheck.200.mc15_8TeV.361106.ZeeJet.pPb_myOutput.root")'&

# ... and for period B
root -l -b -q '../AnalysisCode/RtrkComparison.C (0, 0, true, false, "user.jeouelle.2.4.30hi.calibcheck.200.mc15_8TeV.361107.ZmumuJet.Pbp_myOutput.root")'&
root -l -b -q '../AnalysisCode/RtrkComparison.C (-6, 0, true, false, "user.jeouelle.2.4.30hi.calibcheck.200.mc15_8TeV.361106.ZeeJet.Pbp_myOutput.root")'&

wait

root -l -b -q '../AnalysisCode/RtrkComparisonHist.C ()'
