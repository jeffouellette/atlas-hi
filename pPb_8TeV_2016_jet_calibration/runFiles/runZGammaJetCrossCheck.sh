# 2016 8TeV pPb main physics stream
root -l -b -q '../AnalysisCode/ZGammaJetCrossCheck.C (313063, 0.03, 1, false, true)'&
root -l -b -q '../AnalysisCode/ZGammaJetCrossCheck.C (313067, 1.24, 1, false, true)'&
root -l -b -q '../AnalysisCode/ZGammaJetCrossCheck.C (313100, 9.66, 1, false, true)'&
root -l -b -q '../AnalysisCode/ZGammaJetCrossCheck.C (313107, 11.92, 1, false, true)'& 
root -l -b -q '../AnalysisCode/ZGammaJetCrossCheck.C (313136, 10.4, 1, false, true)'& 
root -l -b -q '../AnalysisCode/ZGammaJetCrossCheck.C (313187, 3.67, 1, false, true)'&

wait

root -l -b -q '../AnalysisCode/ZGammaJetCrossCheck.C (313259, 5.12, 1, false, true)'&
root -l -b -q '../AnalysisCode/ZGammaJetCrossCheck.C (313285, 4.74, 1, false, true)'&
root -l -b -q '../AnalysisCode/ZGammaJetCrossCheck.C (313295, 10.69, 1, false, true)'&
root -l -b -q '../AnalysisCode/ZGammaJetCrossCheck.C (313333, 4.13, 1, false, true)'&
root -l -b -q '../AnalysisCode/ZGammaJetCrossCheck.C (313435, 0.39, 1, false, true)'&
root -l -b -q '../AnalysisCode/ZGammaJetCrossCheck.C (313572, 0.01, 1, false, false)'&

wait

root -l -b -q '../AnalysisCode/ZGammaJetCrossCheck.C (313574, 1.33, 1, false, false)'&
root -l -b -q '../AnalysisCode/ZGammaJetCrossCheck.C (313575, 7.54, 1, false, false)'&
root -l -b -q '../AnalysisCode/ZGammaJetCrossCheck.C (313603, 8.69, 1, false, false)'& 
root -l -b -q '../AnalysisCode/ZGammaJetCrossCheck.C (313629, 6.86, 1, false, false)'&
root -l -b -q '../AnalysisCode/ZGammaJetCrossCheck.C (313630, 7.90, 1, false, false)'& 
root -l -b -q '../AnalysisCode/ZGammaJetCrossCheck.C (313688, 7.96, 1, false, false)'&

wait

root -l -b -q '../AnalysisCode/ZGammaJetCrossCheck.C (313695, 4.53, 1, false, false)'& 
root -l -b -q '../AnalysisCode/ZGammaJetCrossCheck.C (313833, 5.11, 1, false, false)'&
root -l -b -q '../AnalysisCode/ZGammaJetCrossCheck.C (313878, 2.16, 1, false, false)'& 
root -l -b -q '../AnalysisCode/ZGammaJetCrossCheck.C (313929, 0.63, 1, false, false)'& 
root -l -b -q '../AnalysisCode/ZGammaJetCrossCheck.C (313935, 10.96, 1, false, false)'&
root -l -b -q '../AnalysisCode/ZGammaJetCrossCheck.C (313984, 2.40, 1, false, false)'&

wait

root -l -b -q '../AnalysisCode/ZGammaJetCrossCheck.C (314014, 7.36, 1, false, false)'& 
root -l -b -q '../AnalysisCode/ZGammaJetCrossCheck.C (314077, 10.19, 1, false, false)'&
root -l -b -q '../AnalysisCode/ZGammaJetCrossCheck.C (314105, 6.50, 1, false, false)'&
root -l -b -q '../AnalysisCode/ZGammaJetCrossCheck.C (314112, 10.49, 1, false, false)'&
root -l -b -q '../AnalysisCode/ZGammaJetCrossCheck.C (314157, 9.83, 1, false, false)'& 
root -l -b -q '../AnalysisCode/ZGammaJetCrossCheck.C (314170, 4.92, 1, false, false)'& 

wait

# gamma+jet pp sample with pPb data overlay
root -l -b -q '../AnalysisCode/ZGammaJetCrossCheck.C (1, 0, 0.000028838, true, true, "user.jeouelle.2.4.30hi.calibcheck.900.mc15_8TeV.423100.Slice1.pPb_myOutput.root")'&
root -l -b -q '../AnalysisCode/ZGammaJetCrossCheck.C (2, 0, 0.000026548, true, true, "user.jeouelle.2.4.30hi.calibcheck.900.mc15_8TeV.423101.Slice2.pPb_myOutput.root")'&
root -l -b -q '../AnalysisCode/ZGammaJetCrossCheck.C (3, 0, 0.000029816, true, true, "user.jeouelle.2.4.30hi.calibcheck.900.mc15_8TeV.423102.Slice3.pPb_myOutput.root")'&
root -l -b -q '../AnalysisCode/ZGammaJetCrossCheck.C (4, 0, 0.000041955, true, true, "user.jeouelle.2.4.30hi.calibcheck.900.mc15_8TeV.423103.Slice4.pPb_myOutput.root")'&
root -l -b -q '../AnalysisCode/ZGammaJetCrossCheck.C (5, 0, 0.000050105, true, true, "user.jeouelle.2.4.30hi.calibcheck.900.mc15_8TeV.423104.Slice5.pPb_myOutput.root")'&
root -l -b -q '../AnalysisCode/ZGammaJetCrossCheck.C (6, 0, 0.000051776, true, true, "user.jeouelle.2.4.30hi.calibcheck.900.mc15_8TeV.423105.Slice6.pPb_myOutput.root")'&

wait

root -l -b -q '../AnalysisCode/ZGammaJetCrossCheck.C (1, 0, 0.00002885, true, false, "user.jeouelle.2.4.30hi.calibcheck.900.mc15_8TeV.423100.Slice1.Pbp_myOutput.root")'&
root -l -b -q '../AnalysisCode/ZGammaJetCrossCheck.C (2, 0, 0.000026538, true, false, "user.jeouelle.2.4.30hi.calibcheck.900.mc15_8TeV.423101.Slice2.Pbp_myOutput.root")'&
root -l -b -q '../AnalysisCode/ZGammaJetCrossCheck.C (3, 0, 0.00002988, true, false, "user.jeouelle.2.4.30hi.calibcheck.900.mc15_8TeV.423102.Slice3.Pbp_myOutput.root")'&
root -l -b -q '../AnalysisCode/ZGammaJetCrossCheck.C (4, 0, 0.000041994, true, false, "user.jeouelle.2.4.30hi.calibcheck.900.mc15_8TeV.423103.Slice4.Pbp_myOutput.root")'&
root -l -b -q '../AnalysisCode/ZGammaJetCrossCheck.C (5, 0, 0.000050105, true, false, "user.jeouelle.2.4.30hi.calibcheck.900.mc15_8TeV.423104.Slice5.Pbp_myOutput.root")'&
root -l -b -q '../AnalysisCode/ZGammaJetCrossCheck.C (6, 0, 0.000051791, true, false, "user.jeouelle.2.4.30hi.calibcheck.900.mc15_8TeV.423105.Slice6.Pbp_myOutput.root")'&

wait

## validation signal-only gamma+jet sample (period A only!)
#root -l -b -q '../AnalysisCode/ZGammaJetCrossCheck.C (1, 0, true, true, "user.jeouelle.2.4.30hi.calibcheck.valid.420.mc15_8TeV.423100.Slice1.pPb")'&
#root -l -b -q '../AnalysisCode/ZGammaJetCrossCheck.C (2, 0, true, true, "user.jeouelle.2.4.30hi.calibcheck.valid.420.mc15_8TeV.423101.Slice2.pPb")'&
#root -l -b -q '../AnalysisCode/ZGammaJetCrossCheck.C (3, 0, true, true, "user.jeouelle.2.4.30hi.calibcheck.valid.420.mc15_8TeV.423102.Slice3.pPb")'&
#root -l -b -q '../AnalysisCode/ZGammaJetCrossCheck.C (4, 0, true, true, "user.jeouelle.2.4.30hi.calibcheck.valid.420.mc15_8TeV.423103.Slice4.pPb")'&
#root -l -b -q '../AnalysisCode/ZGammaJetCrossCheck.C (5, 0, true, true, "user.jeouelle.2.4.30hi.calibcheck.valid.420.mc15_8TeV.423104.Slice5.pPb")'&
#root -l -b -q '../AnalysisCode/ZGammaJetCrossCheck.C (6, 0, true, true, "user.jeouelle.2.4.30hi.calibcheck.valid.420.mc15_8TeV.423105.Slice6.pPb")'&
#
#wait

# Zmumu (data overlay) and Zee (signal only) samples for period A
root -l -b -q '../AnalysisCode/ZGammaJetCrossCheck.C (0, 0, 1, true, true, "user.jeouelle.2.4.30hi.calibcheck.900.mc15_8TeV.361107.ZmumuJet.pPb_myOutput.root")'&
root -l -b -q '../AnalysisCode/ZGammaJetCrossCheck.C (-6, 0, 1, true, true, "user.jeouelle.2.4.30hi.calibcheck.900.mc15_8TeV.361106.ZeeJet.pPb_myOutput.root")'&
#root -l -b -q '../AnalysisCode/ZGammaJetCrossCheck.C (-1, 0, true, true, "user.jeouelle.2.4.30hi.calibcheck.420.mc15_8TeV.361106.ZeeJet1.pPb_myOutput.root")'&
#root -l -b -q '../AnalysisCode/ZGammaJetCrossCheck.C (-2, 0, true, true, "user.jeouelle.2.4.30hi.calibcheck.420.mc15_8TeV.361106.ZeeJet2.pPb_myOutput.root")'&
#root -l -b -q '../AnalysisCode/ZGammaJetCrossCheck.C (-3, 0, true, true, "user.jeouelle.2.4.30hi.calibcheck.420.mc15_8TeV.361106.ZeeJet3.pPb_myOutput.root")'&
#root -l -b -q '../AnalysisCode/ZGammaJetCrossCheck.C (-4, 0, true, true, "user.jeouelle.2.4.30hi.calibcheck.420.mc15_8TeV.361106.ZeeJet4.pPb_myOutput.root")'&
#root -l -b -q '../AnalysisCode/ZGammaJetCrossCheck.C (-5, 0, true, true, "user.jeouelle.2.4.30hi.calibcheck.420.mc15_8TeV.361106.ZeeJet5.pPb_myOutput.root")'&

wait

# ... and for period B
root -l -b -q '../AnalysisCode/ZGammaJetCrossCheck.C (0, 0, 1, true, false, "user.jeouelle.2.4.30hi.calibcheck.900.mc15_8TeV.361107.ZmumuJet.Pbp_myOutput.root")'&
root -l -b -q '../AnalysisCode/ZGammaJetCrossCheck.C (-6, 0, 1, true, false, "user.jeouelle.2.4.30hi.calibcheck.900.mc15_8TeV.361106.ZeeJet.Pbp_myOutput.root")'&
#root -l -b -q '../AnalysisCode/ZGammaJetCrossCheck.C (-1, 0, true, false, "user.jeouelle.2.4.30hi.calibcheck.420.mc15_8TeV.361106.ZeeJet1.Pbp_myOutput.root")'&
#root -l -b -q '../AnalysisCode/ZGammaJetCrossCheck.C (-2, 0, true, false, "user.jeouelle.2.4.30hi.calibcheck.420.mc15_8TeV.361106.ZeeJet2.Pbp_myOutput.root")'&
#root -l -b -q '../AnalysisCode/ZGammaJetCrossCheck.C (-3, 0, true, false, "user.jeouelle.2.4.30hi.calibcheck.420.mc15_8TeV.361106.ZeeJet3.Pbp_myOutput.root")'&
#root -l -b -q '../AnalysisCode/ZGammaJetCrossCheck.C (-4, 0, true, false, "user.jeouelle.2.4.30hi.calibcheck.420.mc15_8TeV.361106.ZeeJet4.Pbp_myOutput.root")'&
#root -l -b -q '../AnalysisCode/ZGammaJetCrossCheck.C (-5, 0, true, false, "user.jeouelle.2.4.30hi.calibcheck.420.mc15_8TeV.361106.ZeeJet5.Pbp_myOutput.root")'&

wait

root -l -b -q '../AnalysisCode/ZGammaJetCrossCheckHist.C ()'
#root -l -b -q '../AnalysisCode/GammaMCValidOverlayComp.C ()'
