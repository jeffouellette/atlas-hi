version=321

#./bin/runJetTrackEtaPhi xAOD_${version} 313107 false true &
#./bin/runJetTrackEtaPhi xAOD_${version} 313136 false true &
./bin/runJetTrackEtaPhi xAOD_${version} 313295 false true &
./bin/runJetTrackEtaPhi xAOD_${version} 313935 false false &
#./bin/runJetTrackEtaPhi xAOD_${version} 314077 false false &

wait

#version=310
#
#./bin/runJetTrackEtaPhi mc_Ds_${version} 1 true true user.jeouelle.2.4.30hi.jetTreeMaker.310.mc15_8TeV.420185.jetjet.JZ1WA.DstarMinus_minPt1GeV.pPb_myOutput.root 7.307e4 0.00008402 58652 &
#./bin/runJetTrackEtaPhi mc_Ds_${version} 1 true true user.jeouelle.2.4.30hi.jetTreeMaker.310.mc15_8TeV.420186.jetjet.JZRW1B.DstarMinus_mPt1GeV.pPb_myOutput.root 1.8373e4 0.00018431 58576 &
#./bin/runJetTrackEtaPhi mc_Ds_${version} 1 true true user.jeouelle.2.4.30hi.jetTreeMaker.310.mc15_8TeV.420188.jetjet.JZRW1B.DstarPlus_minPt1GeV.pPb_myOutput.root 1.8377e4 0.00018475 98609 &
#
#./bin/runJetTrackEtaPhi mc_Ds_${version} 1 true false user.jeouelle.2.4.30hi.jetTreeMaker.310.mc15_8TeV.420183.jetjet.JZ1WA.DstarPlus_minPt1GeV.Pbp_myOutput.root 7.307e4 0.000084366 458654 &
#./bin/runJetTrackEtaPhi mc_Ds_${version} 1 true false user.jeouelle.2.4.30hi.jetTreeMaker.310.mc15_8TeV.420185.jetjet.JZ1WA.DstarMinus_minPt1GeV.Pbp_myOutput.root 7.307e4 0.000084227 225793 &
#
#wait

./bin/runJetTrackEtaPhiHist
