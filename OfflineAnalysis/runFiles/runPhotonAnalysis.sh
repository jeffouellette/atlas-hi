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
#bin/runPhotonAnalysis 365763 false "Express" &
#bin/runPhotonAnalysis 365768 false "Express" &
#bin/runPhotonAnalysis 365834 false "HardProbes" &
#bin/runPhotonAnalysis 365914 false "HardProbes" &
#
#wait
#
#bin/runPhotonAnalysis 365932 false "HardProbes" &
#bin/runPhotonAnalysis 366011 false "HardProbes" &
#bin/runPhotonAnalysis 366029 false "HardProbes" &
#bin/runPhotonAnalysis 366092 false "HardProbes" &
#bin/runPhotonAnalysis 366142 false "HardProbes" &
#bin/runPhotonAnalysis 366268 false "HardProbes" &
#
#wait
#
#bin/runPhotonAnalysis 366337 false "HardProbes" &
#bin/runPhotonAnalysis 366383 false "HardProbes" &
#bin/runPhotonAnalysis 366413 false "HardProbes" &
#bin/runPhotonAnalysis 366476 false "HardProbes" &
#bin/runPhotonAnalysis 366526 false "HardProbes" &
#bin/runPhotonAnalysis 366528 false "HardProbes" &
#
#wait
#
#bin/runPhotonAnalysis 366627 false "HardProbes" &
#bin/runPhotonAnalysis 366691 false "HardProbes" &
#bin/runPhotonAnalysis 366754 false "HardProbes" &
#bin/runPhotonAnalysis 366805 false "HardProbes" &
#bin/runPhotonAnalysis 366860 false "HardProbes" &
#bin/runPhotonAnalysis 366878 false "HardProbes" &
#
#wait
#
#bin/runPhotonAnalysis 366919 false "HardProbes" &
#bin/runPhotonAnalysis 366931 false "HardProbes" &
#bin/runPhotonAnalysis 366994 false "HardProbes" &
#bin/runPhotonAnalysis 367023 false "HardProbes" &
#bin/runPhotonAnalysis 367099 false "HardProbes" &
bin/runPhotonAnalysis 367134 false "HardProbes" &

wait

bin/runPhotonAnalysis 367165 false "HardProbes" &
bin/runPhotonAnalysis 367170 false "HardProbes" &
#bin/runPhotonAnalysis 367233 false "Express" &
bin/runPhotonAnalysis 367273 false "HardProbes" &
#bin/runPhotonAnalysis 367318 false "Express" &
bin/runPhotonAnalysis 367321 false "HardProbes" &

wait

bin/runPhotonAnalysis 367363 false "HardProbes" &
bin/runPhotonAnalysis 367364 false "HardProbes" &
bin/runPhotonAnalysis 367365 false "HardProbes" &
bin/runPhotonAnalysis 367384 false "HardProbes" &

wait

rm ../RootFiles/PhotonAnalysis/outFile.root
hadd ../RootFiles/PhotonAnalysis/outFile.root ../RootFiles/PhotonAnalysis/dataSet_36*.root

bin/runPhotonAnalysisHist
