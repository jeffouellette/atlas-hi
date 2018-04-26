#include "../GlobalParams.C"

/** User defined parameters **/

const int useDataVersion = 8; // Specifies which version of raw data to use. Different versions have different branches in trees and (generally) different analysis procedures. If any value is given other than the currently accepted value, every run will be skipped and nothing will happen. If 0 is specified, MC will be run instead. Currently accepted values: 0, 8.
const bool runPeriodA = true; // Analyze period A data
const bool runPeriodB = true; // Analyze period B data
const bool scaleAnalyses = true; // Whether to separate analysis plots by relatively rescaling each plot
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

const double pbins[68] = {15., 20., 25., 30., 35., 40., 45., 50., 55., 60., 65., 70., 75., 80., 85., 90., 95., 100., 105., 110., 115., 120., 125., 130., 135., 140., 145., 150., 155., 160., 165., 170., 175., 180., 185., 190., 195., 200., 205., 210., 220., 230., 240., 250., 260., 270., 280., 290., 300., 320., 340., 360., 380., 400., 425., 450., 475., 500., 550., 600., 700., 800., 1000., 1250., 1500., 2000., 2500., 6000.};
const int numpbins = sizeof(pbins)/sizeof(pbins[0]) - 1;
const double etabins[7] = {-1.3, -0.8, -0.3, 0, 0.3, 0.8, 1.3};
//const double etabins[2] = {-4.9, 4.9}; // Used for avoiding eta-binning.
const int numetabins = sizeof(etabins)/sizeof(etabins[0]) - 1;

/** End general parameters **/
