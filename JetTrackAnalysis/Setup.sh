#! /bin/bash

export LD_LIBRARY_PATH=$ATLAS_PATH/AnalysisCode/lib:$JETTRACKS_PATH/AnalysisCode/lib:$LD_LIBRARY_PATH
export DYLD_LIBRARY_PATH=$ATLAS_PATH/AnalysisCode/lib:$JETTRACKS_PATH/AnalysisCode/lib:$DYLD_LIBRARY_PATH

alias g++=/usr/bin/clang++
alias gcc=/usr/bin/clang
