#include "../GlobalParams.C"

/** User defined parameters **/

const int useDataVersion = 8; // Specifies which version of raw data to use. Different versions have different branches in trees and (generally) different analysis procedures. If any value is given other than the currently accepted value, every run will be skipped and nothing will happen. If 0 is specified, MC will be run instead. Currently accepted values: 0, 8.
const bool runPeriodA = true; // Analyze period A data
const bool runPeriodB = true; // Analyze period B data
const bool scaleAnalyses = true; // Whether to separate analysis plots by relatively rescaling each plot
workPath = homePath + "P_Pb_dijet_analyses/"; // Home analysis directory, should be modified in code outside this path structure
const string fittedFunctionType = /*"erf";*/ "fermi_dirac"; // Function used to fit trigger efficiencies- accepted values: "fermi_dirac", "erf"

const double dijetPtRatioCut = 0.4; // Maximum pt allowed for subsubleading jet as a proportion of the leading jet pt
const double dijetMinimumPt = 20; // Minimum Pt required for jets in dijet analysis
const double thirdJetMaximumPt = 10; // Maximum Pt allowed for subsubleading jet
const bool highPtJetsOnly = false; // only use jets above 70GeV

//const int run_list_v3[30] = {313063, 313067, 313100, 313107, 313136, 313187, 313259, 313285, 313295, 313333, 313435, 313572, 313574, 313575, 313603, 313629, 313630, 313688, 313695, 313833, 313878, 313929, 313935, 313984, 314014, 314077, 314105, 314112, 314157, 314170}; // full run list for future reference
const int run_list_v8[30] = {313063, 313067, 313100, 313107, 313136, 313187, 313259, 313285, 313295, 313333, 313435, 313572, 313574, 313575, 313603, 313629, 313630, 313688, 313695, 313833, 313878, 313929, 313935, 313984, 314014, 314077, 314105, 314112, 314157, 314170};
const int mc_samples[9] = {13202455, 13211348, 13211339, 13216170, 13214408, 13211331, 13200714, 13237282, 13237282}; // these are the task id's attached to each MC sample

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
