#version=313
#
#./bin/runJetTrackEtaPhi xAOD_${version} 313107 false true &
#./bin/runJetTrackEtaPhi xAOD_${version} 313136 false true &
#./bin/runJetTrackEtaPhi xAOD_${version} 313295 false true &
#
#./bin/runJetTrackEtaPhi xAOD_${version} 313935 false false &
#./bin/runJetTrackEtaPhi xAOD_${version} 314077 false false &
##./bin/runJetTrackEtaPhi xAOD_${version} 314112 false false &
#
#wait

version=310

#./bin/runJetTrackEtaPhi mc_Ds_${version} 1 true true user.jeouelle.2.4.30hi.jetTreeMaker.310.mc15_8TeV.420183.jetjet.JZ1WA.DstarPlus_minPt1GeV.pPb_myOutput.root 7.307e4 0.00008435 95576 &  #mc15_pPb8TeV.420183.Pythia8EvtGen_A14NNPDF23LO_jetjet_JZ1WA_DstarPlus_minPt1GeV.merge.AOD.e6181_e5984_d1439_r9645_r9647
./bin/runJetTrackEtaPhi mc_Ds_${version} 1 true true user.jeouelle.2.4.30hi.jetTreeMaker.310.mc15_8TeV.420185.jetjet.JZ1WA.DstarMinus_minPt1GeV.pPb_myOutput.root 7.307e4 0.00008402 58652 &  #mc15_pPb8TeV.420185.Pythia8EvtGen_A14NNPDF23LO_jetjet_JZ1WA_DstarMinus_minPt1GeV.merge.AOD.e6181_e5984_d1439_r9645_r9647
./bin/runJetTrackEtaPhi mc_Ds_${version} 1 true true user.jeouelle.2.4.30hi.jetTreeMaker.310.mc15_8TeV.420186.jetjet.JZRW1B.DstarMinus_minPt1GeV.pPb_myOutput.root 1.8373e4 0.00018431 58576 &  #mc15_pPb8TeV.420186.Pythia8EvtGen_A14NNPDF23LO_jetjet_JZRW1B_DstarMinus_mPt1GeV.merge.AOD.e6181_e5984_d1439_r9645_r9647
./bin/runJetTrackEtaPhi mc_Ds_${version} 1 true true user.jeouelle.2.4.30hi.jetTreeMaker.310.mc15_8TeV.420188.jetjet.JZRW1B.DstarPlus_minPt1GeV.pPb_myOutput.root 1.8377e4 0.00018475 98609 &  #mc15_pPb8TeV.420188.Pythia8EvtGen_A14NNPDF23LO_jetjet_JZRW1B_DstarPlus_minPt1GeV.merge.AOD.e6181_e5984_d1439_r9645_r9647

./bin/runJetTrackEtaPhi mc_Ds_${version} 1 true false user.jeouelle.2.4.30hi.jetTreeMaker.310.mc15_8TeV.420183.jetjet.JZ1WA.DstarPlus_minPt1GeV.Pbp_myOutput.root 7.307e4 0.000084366 458654 &  #mc15_pPb8TeV.420183.Pythia8EvtGen_A14NNPDF23LO_jetjet_JZ1WA_DstarPlus_minPt1GeV.merge.AOD.e6182_e5984_d1440_r9645_r9647
./bin/runJetTrackEtaPhi mc_Ds_${version} 1 true false user.jeouelle.2.4.30hi.jetTreeMaker.310.mc15_8TeV.420185.jetjet.JZ1WA.DstarMinus_minPt1GeV.Pbp_myOutput.root 7.307e4 0.000084227 225793 &  #mc15_pPb8TeV.420185.Pythia8EvtGen_A14NNPDF23LO_jetjet_JZ1WA_DstarMinus_minPt1GeV.merge.AOD.e6182_e5984_d1440_r9645_r9647
./bin/runJetTrackEtaPhi mc_Ds_${version} 1 true false user.jeouelle.2.4.30hi.jetTreeMaker.310.mc15_8TeV.420186.jetjet.JZRW1B.DstarMinus_minPt1GeV.Pbp_myOutput.root 1.8374e4 0.00018322 225793 &  #mc15_pPb8TeV.420186.Pythia8EvtGen_A14NNPDF23LO_jetjet_JZRW1B_DstarMinus_mPt1GeV.merge.AOD.e6182_e5984_d1440_r9645_r9647
./bin/runJetTrackEtaPhi mc_Ds_${version} 1 true false user.jeouelle.2.4.30hi.jetTreeMaker.310.mc15_8TeV.420188.jetjet.JZRW1B.DstarPlus_minPt1GeV.pPb_myOutput.root 1.8377e4 0.00018475 98609 &  #mc15_pPb8TeV.420188.Pythia8EvtGen_A14NNPDF23LO_jetjet_JZRW1B_DstarPlus_minPt1GeV.merge.AOD.e6182_e5984_d1440_r9645_r9647

wait

./bin/runJetTrackEtaPhiHist
