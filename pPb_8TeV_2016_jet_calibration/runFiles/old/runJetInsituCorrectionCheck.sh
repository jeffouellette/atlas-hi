# 2016 8TeV pPb main physics stream
root -l -b -q '../AnalysisCode/JetInsituCorrectionCheck.C (313063, 0.03, true)'&
root -l -b -q '../AnalysisCode/JetInsituCorrectionCheck.C (313067, 1.24, true)'&
root -l -b -q '../AnalysisCode/JetInsituCorrectionCheck.C (313100, 9.66, true)'&
root -l -b -q '../AnalysisCode/JetInsituCorrectionCheck.C (313107, 11.92, true)'& 
root -l -b -q '../AnalysisCode/JetInsituCorrectionCheck.C (313136, 10.4, true)'& 
root -l -b -q '../AnalysisCode/JetInsituCorrectionCheck.C (313187, 3.67, true)'&

wait

root -l -b -q '../AnalysisCode/JetInsituCorrectionCheck.C (313259, 5.12, true)'&
root -l -b -q '../AnalysisCode/JetInsituCorrectionCheck.C (313285, 4.74, true)'&
root -l -b -q '../AnalysisCode/JetInsituCorrectionCheck.C (313295, 10.69, true)'&
root -l -b -q '../AnalysisCode/JetInsituCorrectionCheck.C (313333, 4.13, true)'&
root -l -b -q '../AnalysisCode/JetInsituCorrectionCheck.C (313435, 0.39, true)'&
root -l -b -q '../AnalysisCode/JetInsituCorrectionCheck.C (313572, 0.01, false)'&

wait

root -l -b -q '../AnalysisCode/JetInsituCorrectionCheck.C (313574, 1.33, false)'&
root -l -b -q '../AnalysisCode/JetInsituCorrectionCheck.C (313575, 7.54, false)'&
root -l -b -q '../AnalysisCode/JetInsituCorrectionCheck.C (313603, 8.69, false)'& 
root -l -b -q '../AnalysisCode/JetInsituCorrectionCheck.C (313629, 6.86, false)'&
root -l -b -q '../AnalysisCode/JetInsituCorrectionCheck.C (313630, 7.90, false)'& 
root -l -b -q '../AnalysisCode/JetInsituCorrectionCheck.C (313688, 7.96, false)'&

wait

root -l -b -q '../AnalysisCode/JetInsituCorrectionCheck.C (313695, 4.53, false)'& 
root -l -b -q '../AnalysisCode/JetInsituCorrectionCheck.C (313833, 5.11, false)'&
root -l -b -q '../AnalysisCode/JetInsituCorrectionCheck.C (313878, 2.16, false)'& 
root -l -b -q '../AnalysisCode/JetInsituCorrectionCheck.C (313929, 0.63, false)'& 
root -l -b -q '../AnalysisCode/JetInsituCorrectionCheck.C (313935, 10.96, false)'&
root -l -b -q '../AnalysisCode/JetInsituCorrectionCheck.C (313984, 2.40, false)'&

wait

root -l -b -q '../AnalysisCode/JetInsituCorrectionCheck.C (314014, 7.36, false)'& 
root -l -b -q '../AnalysisCode/JetInsituCorrectionCheck.C (314077, 10.19, false)'&
root -l -b -q '../AnalysisCode/JetInsituCorrectionCheck.C (314105, 6.50, false)'&
root -l -b -q '../AnalysisCode/JetInsituCorrectionCheck.C (314112, 10.49, false)'&
root -l -b -q '../AnalysisCode/JetInsituCorrectionCheck.C (314157, 9.83, false)'& 
root -l -b -q '../AnalysisCode/JetInsituCorrectionCheck.C (314170, 4.92, false)'& 

wait

root -l -b -q '../AnalysisCode/JetInsituCorrectionCheckHist.C ()'
