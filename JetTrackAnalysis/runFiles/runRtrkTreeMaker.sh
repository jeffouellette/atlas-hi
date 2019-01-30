#version=009
#
#./bin/runRtrkTreeMaker mc_DijetSignal_${version} 1 true true user.jeouelle.2.4.30hi.jetTreeMaker.signalonly.${version}.mc15_8TeV.420011.jetjet.JZ1R04.pPb_myOutput.root 1 1 1 & # mc15_pPb8TeV.420011.Pythia8EvtGen_A14NNPDF23LO_jetjet_JZ1R04.merge.AOD.e6518_s3084_s3153_r9985_r9647
#./bin/runRtrkTreeMaker mc_DijetSignal_${version} 2 true true user.jeouelle.2.4.30hi.jetTreeMaker.signalonly.${version}.mc15_8TeV.420012.jetjet.JZ2R04.pPb_myOutput.root 1 1 1 & # mc15_pPb8TeV.420012.Pythia8EvtGen_A14NNPDF23LO_jetjet_JZ2R04.merge.AOD.e6394_s3084_s3153_r9985_r9647
##./bin/runRtrkTreeMaker mc_DijetSignal_${version} 3 true true user.jeouelle.2.4.30hi.jetTreeMaker.signalonly.${version}.mc15_8TeV.420013.jetjet.JZ3R04.pPb_myOutput.root 1 1 1 & # mc15_pPb8TeV.420013.Pythia8EvtGen_A14NNPDF23LO_jetjet_JZ3R04.merge.AOD.e6518_s3084_s3153_r9985_r9647
##./bin/runRtrkTreeMaker mc_DijetSignal_${version} 4 true true user.jeouelle.2.4.30hi.jetTreeMaker.signalonly.${version}.mc15_8TeV.420014.jetjet.JZ4R04.pPb_myOutput.root 1 1 1 & # mc15_pPb8TeV.420014.Pythia8EvtGen_A14NNPDF23LO_jetjet_JZ4R04.merge.AOD.e6518_s3084_s3153_r9985_r9647
##./bin/runRtrkTreeMaker mc_DijetSignal_${version} 5 true true user.jeouelle.2.4.30hi.jetTreeMaker.signalonly.${version}.mc15_8TeV.420015.jetjet.JZ5R04.pPb_myOutput.root 1 1 1 & # mc15_pPb8TeV.420015.Pythia8EvtGen_A14NNPDF23LO_jetjet_JZ5R04.merge.AOD.e6518_s3084_s3153_r9985_r9647

wait

version=011
#./bin/runRtrkTreeMaker mc_DijetOverlay_${version} 1 true true user.jeouelle.2.4.30hi.jetTreeMaker.${version}.mc15_8TeV.420011.jetjet.JZ1R04.pPb_myOutput.root 1 1 1 & # mc15_pPb8TeV.420011.Pythia8EvtGen_A14NNPDF23LO_jetjet_JZ1R04.merge.AOD.e6518_d1489_r10779_r9647
./bin/runRtrkTreeMaker mc_DijetOverlay_${version} 2 true true user.jeouelle.2.4.30hi.jetTreeMaker.${version}.mc15_8TeV.420012.jetjet.JZ2R04.pPb_myOutput.root 1 1 1 & # mc15_pPb8TeV.420011.Pythia8EvtGen_A14NNPDF23LO_jetjet_JZ1R04.merge.AOD.e6518_d1489_r10779_r9647
#./bin/runRtrkTreeMaker mc_DijetOverlay_${version} 3 true true user.jeouelle.2.4.30hi.jetTreeMaker.${version}.mc15_8TeV.420013.jetjet.JZ3R04.pPb_myOutput.root 1 1 1 & # mc15_pPb8TeV.420011.Pythia8EvtGen_A14NNPDF23LO_jetjet_JZ1R04.merge.AOD.e6518_d1489_r10779_r9647
#./bin/runRtrkTreeMaker mc_DijetOverlay_${version} 4 true true user.jeouelle.2.4.30hi.jetTreeMaker.${version}.mc15_8TeV.420014.jetjet.JZ4R04.pPb_myOutput.root 1 1 1 & # mc15_pPb8TeV.420011.Pythia8EvtGen_A14NNPDF23LO_jetjet_JZ1R04.merge.AOD.e6518_d1489_r10779_r9647
#./bin/runJetTrackEtaPhi mc_DijetOverlay_${version} 5 true true user.jeouelle.2.4.30hi.jetTreeMaker.${version}.mc15_8TeV.420015.jetjet.JZ5R04.pPb_myOutput.root 1 1 1 & # mc15_pPb8TeV.420011.Pythia8EvtGen_A14NNPDF23LO_jetjet_JZ1R04.merge.AOD.e6518_d1489_r10779_r9647

wait

echo "All done!"
