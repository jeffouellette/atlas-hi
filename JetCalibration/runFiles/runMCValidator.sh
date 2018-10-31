version=430
./bin/runMCValidator jetjet_valid_${version} 2 true 60 160 user.jeouelle.2.4.30hi.calibcheck.valid.${version}.mc15_8TeV.420012.jetjet.JZ2R04.pPb_myOutput.root 1e3 1 9875 &
./bin/runMCValidator jetjet_valid_${version} 2 false 60 160 user.jeouelle.2.4.30hi.calibcheck.valid.${version}.mc15_8TeV.420012.jetjet.JZ2R04.Pbp_myOutput.root 1e3 1 9835 &

wait

./bin/runMCValidatorHist
