#! /bin/bash

export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$ATLAS_PATH/AnalysisCode/lib:$JETTRACKS_PATH/AnalysisCode/lib
export DYLD_LIBRARY_PATH=$DYLD_LIBRARY_PATH:$ATLAS_PATH/AnalysisCode/lib:$JETTRACKS_PATH/AnalysisCode/lib

alias g++=/usr/bin/clang++
alias gcc=/usr/bin/clang
