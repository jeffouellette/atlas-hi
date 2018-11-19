bin/runPhotonAnalysis 365498 false 0.2674 "HardProbes" &
bin/runPhotonAnalysis 365502 false 1.018 "HardProbes" &
bin/runPhotonAnalysis 365512 false 0.9223 "HardProbes" &
bin/runPhotonAnalysis 365573 false 7.852 "HardProbes" &
bin/runPhotonAnalysis 365602 false 9.897 "HardProbes" &
bin/runPhotonAnalysis 365627 false 25.151 "HardProbes" &

wait

bin/runPhotonAnalysis 365678 false 16.219 "HardProbes" &
bin/runPhotonAnalysis 365681 false 24.847 "HardProbes" &
#bin/runPhotonAnalysis 365709 false 35.472 "HardProbes" &
#bin/runPhotonAnalysis 365752 false 36.751 "HardProbes" &
#bin/runPhotonAnalysis 365763 false 5.484 "HardProbes" &
#bin/runPhotonAnalysis 365768 false 7.659 "HardProbes" &

#wait

bin/runPhotonAnalysis 365834 false 5.651 "HardProbes" &
#bin/runPhotonAnalysis 365914 false 6.595 "HardProbes" &
#bin/runPhotonAnalysis 365932 false 6.595 "HardProbes" &
#bin/runPhotonAnalysis 366011 false 43.806 "HardProbes" &
#bin/runPhotonAnalysis 366029 false 42.881 "HardProbes" &
#bin/runPhotonAnalysis 366092 false 45.336 "HardProbes" &

wait

#bin/runPhotonAnalysis 366142 false 41.946 "HardProbes" &
#bin/runPhotonAnalysis 366268 false 54.532 "HardProbes" &
#bin/runPhotonAnalysis 366337 false  "HardProbes" &

rm ../RootFiles/PhotonAnalysis/outFile.root
hadd ../RootFiles/PhotonAnalysis/outFile.root ../RootFiles/PhotonAnalysis/dataSet_36*.root

bin/runPhotonAnalysisHist
