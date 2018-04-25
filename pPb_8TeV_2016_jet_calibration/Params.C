#include "../GlobalParams.C"

/** User defined parameters **/

const bool runPeriodA = true; // Analyze period A data
const bool runPeriodB = true; // Analyze period B data
const bool scaleAnalyses = true; // Whether to separate analysis plots by relatively rescaling each plot
workPath = homePath + "pPb_8TeV_2016_jet_calibration/"; // Home analysis directory, should be modified in code outside this path structure

/** End user defined parameters **/


/** General (non-user defined) paramters **/

// More directory information - PLEASE DO NOT CHANGE!!! These values are overwritten when calling triggerUtil::initialize().
rootPath = workPath + "rootFiles/"; // Where analyzed *.root files are stored. Different analysis modules have different subdirectories here.
dataPath = workPath + "data/"; // Where the *.root raw data files (from the CERN grid) are stored.
plotPath = workPath + "Plots/"; // Where plots are stored.
ptPath = rootPath + "ptData/"; // Where the pt analysis module output is stored.
trigPath = rootPath + "trigData/";  // Where the trigger fire count module output is stored.
effPath = rootPath + "effData/"; // Where the trigger efficiency module output is stored.
xPath = rootPath + "xData/"; // Where the xa/xp module output is stored.
RpPbPath = rootPath + "RpPbData/"; // Where the R_pPb module output is stored.

/** End general parameters **/
