/**
 * Contains useful variables and directory information for 2016 pPb data analyses.
 * Author: Jeff Ouellette
 * Dated: 4/25/2018
 */

/** User defined parameters **/

static const bool debugStatements = false; // Print out periodic statements to monitor code flow
const bool scaleAnalyses = true; // Whether to separate analysis plots by relatively rescaling each plot
static const string homePath = "/Users/jeffouellette/Research/atlas-hi/"; // ATLAS Heavy Ions home directory

const double dR_HEC = 0.4; // details on the hadronic end cap data cuts.
const double lowerPhiCut = TMath::Pi()-dR_HEC;
const double upperPhiCut = 3.*TMath::Pi()/2.+dR_HEC;
const double lowerEtaCut = 1.5-dR_HEC;
const double upperEtaCut = 3.2+dR_HEC;

const int full_run_list[30] = {313063, 313067, 313100, 313107, 313136, 313187, 313259, 313285, 313295, 313333, 313435, 313572, 313574, 313575, 313603, 313629, 313630, 313688, 313695, 313833, 313878, 313929, 313935, 313984, 314014, 314077, 314105, 314112, 314157, 314170};

/** End user defined parameters **/


/** General (non-user defined) paramters **/

// More directory information - PLEASE DO NOT CHANGE!!! These values are overwritten when calling triggerUtil::initialize().
string workPath; // Home analysis directory, should be modified in code outside this path structure
string rootPath; // Where analyzed *.root files are stored. Different analysis modules have different subdirectories here.
string dataPath; // Where the *.root raw data files (from the CERN grid) are stored.
string plotPath; // Where plots are stored.
string ptPath; // Where the pt analysis module output is stored.
string trigPath;  // Where the trigger fire count module output is stored.
string effPath; // Where the trigger efficiency module output is stored.
string xPath; // Where the xa/xp module output is stored.
string RpPbPath; // Where the R_pPb module output is stored.

// Transverse momentum and pseudorapidity binning
static const int MAX_PT = 6000; // Maximum transverse momentum
static const double MIN_ETA = -4.9; // Minimum detectable pseudorapidity in hadronic calorimeter
static const double MAX_ETA = 4.9; // Maximum detectable pseudorapidity

int numtrigs = 0; // Total number of triggers

// Useful constants
static const float Z = 82;   // value of Z for Pb
static const float A = 208;  // value of A for Pb
static const float sqrt_s_nn = 8160; // Collision energy in CoM frame (GeV)
static const float muon_mass = 0.105658; // GeV

const double pi = TMath::Pi();

/** End general parameters **/


/**
 * Modifies the directory strings to point to the correct locations.
 */
void setupDirectories (const string dataSubDir, const string thisWorkPath) {

  if (thisWorkPath.at(thisWorkPath.length()-1) != '/') workPath = homePath + thisWorkPath + "/";
  else workPath = homePath + thisWorkPath;

  rootPath = workPath + "rootFiles/" + dataSubDir;
  dataPath = workPath + "data/" + dataSubDir;
  plotPath = workPath + "Plots/" + dataSubDir;
  ptPath = rootPath + "ptData/";
  trigPath = rootPath + "trigData/";
  effPath = rootPath + "effData/";
  xPath = rootPath + "xData/";
  RpPbPath = rootPath + "RpPbData/";

  return;
}
