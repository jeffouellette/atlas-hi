#! /bin/bash

export LD_LIBRARY_PATH=$ATLAS_PATH/AnalysisCode/lib:$ROOT_UTIL_PATH/
export DYLD_LIBRARY_PATH=$ATLAS_PATH/AnalysisCode/lib:$ROOT_UTIL_PATH/

alias g++=/usr/bin/clang++
alias gcc=/usr/bin/clang
