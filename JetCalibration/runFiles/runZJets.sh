rm ../rootFiles/ZJets/outFile.root
hadd ../rootFiles/ZJets/outFile.root ../rootFiles/ZeeJets/outFile.root ../rootFiles/ZmumuJets/outFile.root

./bin/runZJetsHist
