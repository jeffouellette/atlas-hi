# atlas-hi #
A git repository for my ATLAS heavy ions analyses
Jeff Ouellette
Last updated: 12/1/2017

## pPb 8TeV 2017 Dijet Analyses ##
Subdirectory allocated towards dijet analyses in the 2016 8.16TeV proton-lead run. This includes two main subdirectories: jet analyses, which runs over the data given a trigger configuration, and trigger analyses, which runs over the triggers and creates a configuration.

## pp 5TeV 2017 Offline Electron Analyses ##
This subdirectory is allocated for offline electron performance analyses in the 2017 5.02TeV proton-proton run. At this point, the subproject is largely finished due to the end of the run. However the code found here may be useful in future offline performance studies.

## pPb 8TeV 2016 Jet Calibration ##
This subdirectory contains code and output from my ATLAS qualification task, which is primarily to perform a jet calibration on pPb 8TeV MC15 data.

## pPb 8TeV 2016 Jet Track Analysis ##
This subdirectory is for analyzing jet-track correlations in the 2016 p+Pb 8.16TeV dataset.

## Using shared libraries in a ROOT interactive session ##
To compile these C++ files, simply run any of the relevant makefiles. This will produce shared libraries in the lib subdirectory. Due to the include interdependencies, the C++ source files cannot be directly run with ROOT as macros. However using ROOT 6 or higher it is possible to load from the shared library directly. To do this, open root:
$  root
root [0] R__LOAD_LIBRARY (/Users/jeffouellette/Research/atlas-hi/lib/libUtils.so)
root [1] R__ADD_INCLUDE_PATH (/Users/jeffouellette/Research/atlas-hi/include)
root [2] #include "Utils.h"
