bin/runPhotonAnalysis 365498 false &
bin/runPhotonAnalysis 365502 false &
bin/runPhotonAnalysis 365512 false &
bin/runPhotonAnalysis 365573 false &
bin/runPhotonAnalysis 365602 false &

wait

rm ../RootFiles/PhotonAnalysis/outFile.root
hadd ../RootFiles/PhotonAnalysis/outFile.root ../RootFiles/PhotonAnalysis/dataSet_36*.root

bin/runPhotonAnalysisHist
