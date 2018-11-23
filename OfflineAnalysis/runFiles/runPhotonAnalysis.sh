#bin/runPhotonAnalysis 365498 false "HardProbes" &
#bin/runPhotonAnalysis 365502 false "HardProbes" &
#bin/runPhotonAnalysis 365512 false "HardProbes" &
#bin/runPhotonAnalysis 365573 false "HardProbes" &
#bin/runPhotonAnalysis 365602 false "HardProbes" &
#bin/runPhotonAnalysis 365627 false "HardProbes" &
#
#wait
#
#bin/runPhotonAnalysis 365678 false "HardProbes" &
#bin/runPhotonAnalysis 365681 false "HardProbes" &
#bin/runPhotonAnalysis 365709 false "HardProbes" &
#bin/runPhotonAnalysis 365752 false "HardProbes" &
#bin/runPhotonAnalysis 365834 false "HardProbes" &
#bin/runPhotonAnalysis 365914 false "HardProbes" &
#
#wait
#
#bin/runPhotonAnalysis 365932 false "HardProbes" &
#bin/runPhotonAnalysis 366011 false "HardProbes" &
#bin/runPhotonAnalysis 366029 false "HardProbes" &
bin/runPhotonAnalysis 366092 false "HardProbes" &
bin/runPhotonAnalysis 366142 false "HardProbes" &
#bin/runPhotonAnalysis 366268 false "HardProbes" &

wait

#bin/runPhotonAnalysis 366337 false "HardProbes" &

rm ../RootFiles/PhotonAnalysis/outFile.root
hadd ../RootFiles/PhotonAnalysis/outFile.root ../RootFiles/PhotonAnalysis/dataSet_36*.root

bin/runPhotonAnalysisHist
