/** User defined parameters **/

const int useDataVersion = 7; // Specifies which version of raw data to use. Different versions have different branches in trees and (generally) different analysis procedures. If any value is given other than the currently accepted value, every run will be skipped and nothing will happen. Currently accepted values: 7
const bool runPeriodA = true; // Analyze period A data
const bool runPeriodB = true; // Analyze period B data
const bool debugStatements = false; // Print out periodic statements to monitor code flow
const bool scaleAnalyses = true; // Whether to separate analysis plots by relatively rescaling each plot
const int trigthres = 0; // Additional threshold requirement for triggers
const string workPath = "/Users/jeffouellette/Research/atlas-hi/P_Pb_dijet_analyses/"; // Home analysis directory
const string minbiasTriggerName = "HLT_mb_mbts_L1MBTS_1"; // Valid triggers: HLT_mb_sptrk_L1MBTS_1, HLT_mb_mbts_L1MBTS_1
const string fittedFunctionType = /*"erf";*/ "fermi_dirac"; // Function used to fit trigger efficiencies- accepted values: "fermi_dirac", "erf"

const double lowerPhiCut = TMath::Pi()-0.4;
const double upperPhiCut = 3.*TMath::Pi()/2.+0.4;
const double lowerEtaCut = 1.5-0.4;
const double upperEtaCut = 3.2+0.4;
const double dijetPtRatioCut = 0.4; // Maximum pt allowed for subsubleading jet as a proportion of the leading jet pt
const double dijetMinimumPt = 20; // Minimum Pt required for jets in dijet analysis
const double thirdJetMaximumPt = 10;

//const int run_list_v3[30] = {313063, 313067, 313100, 313107, 313136, 313187, 313259, 313285, 313295, 313333, 313435, 313572, 313574, 313575, 313603, 313629, 313630, 313688, 313695, 313833, 313878, 313929, 313935, 313984, 314014, 314077, 314105, 314112, 314157, 314170}; // full run list for future reference
const int run_list_v7[26] = {313063, 313067, 313100, 313107, 313136, 313187, 313259, 313285, 313295, 313333, 313435, 313574, 313575, 313629, 313630, 313688, 313695, 313929, 313935, 313984, 314014, 314077, 314105, 314112, 314157, 314170};

/** End user defined parameters **/


/** General (non-user defined) paramters **/

// More directory information - PLEASE DO NOT CHANGE!!! These values are overwritten when calling triggerUtil::initialize().
string rootPath = workPath + "rootFiles/"; // Where analyzed *.root files are stored. Different analysis modules have different subdirectories here.
string dataPath = workPath + "data/"; // Where the *.root raw data files (from the CERN grid) are stored.
string plotPath = workPath + "Plots/"; // Where plots are stored.
string ptPath = rootPath + "ptData/"; // Where the pt analysis module output is stored.
string trigPath = rootPath + "trigData/";  // Where the trigger fire count module output is stored.
string effPath = rootPath + "effData/"; // Where the trigger efficiency module output is stored.
string xPath = rootPath + "xData/"; // Where the xa/xp module output is stored.
string RpPbPath = rootPath + "RpPbData/"; // Where the R_pPb module output is stored.

// Transverse momentum and pseudorapidity binning
const int MAX_PT = 6000; // Maximum transverse momentum
const double MIN_ETA = -4.9; // Minimum detectable pseudorapidity in hadronic calorimeter
const double MAX_ETA = 4.9; // Maximum detectable pseudorapidity

const double pbins[68] = {15., 20., 25., 30., 35., 40., 45., 50., 55., 60., 65., 70., 75., 80., 85., 90., 95., 100., 105., 110., 115., 120., 125., 130., 135., 140., 145., 150., 155., 160., 165., 170., 175., 180., 185., 190., 195., 200., 205., 210., 220., 230., 240., 250., 260., 270., 280., 290., 300., 320., 340., 360., 380., 400., 425., 450., 475., 500., 550., 600., 700., 800., 1000., 1250., 1500., 2000., 2500., 6000.};
const int numpbins = sizeof(pbins)/sizeof(pbins[0]) - 1;
const double etabins[9] = {-4.9, -3.2, -2., -1., 0, 1, 2., 3.2, 4.9};
//const double etabins[2] = {-4.9, 4.9}; // Used for avoiding eta-binning.
const int numetabins = sizeof(etabins)/sizeof(etabins[0]) - 1;
int numtrigs; // Total number of triggers

// Useful constants
const float Z = 82;   // value of Z for Pb
const float A = 208;  // value of A for Pb
const float sqrt_s_nn = 8160; // Collision energy in CoM frame (GeV)

const double pi = TMath::Pi();

/** End general parameters **/
