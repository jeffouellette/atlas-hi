/**
 * Contains useful variables and directory information for 2016 pPb data analyses.
 * Author: Jeff Ouellette
 * Dated: 4/25/2018
 */

/** User defined parameters **/

static const bool debugStatements = false; // Print out periodic statements to monitor code flow
static const string homePath = "/Users/jeffouellette/Research/atlas-hi/"; // ATLAS Heavy Ions home directory
static const string drivePath = "/Volumes/My Passport/Research/atlas-hi/"; // ATLAS Heavy Ions external drive directory

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
string externalWorkPath; // External drive storage directory, should be modified in code below
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
static const double Z = 82;   // value of Z for Pb
static const double A = 208;  // value of A for Pb
static const double sqrt_s_nn = 8160; // Collision energy in CoM frame (GeV)
static const double electron_mass = 0.000511; // GeV
static const double muon_mass = 0.105658; // GeV
static const double Z_mass = 91.2; // GeV

static const double pi = TMath::Pi();

/** End general parameters **/


/**
 * Modifies the directory strings to point to the correct locations.
 */
void setupDirectories (const string dataSubDir, const string thisWorkPath) {

  if (thisWorkPath.at(thisWorkPath.length()-1) != '/') {
    workPath = homePath + thisWorkPath + "/";
    externalWorkPath = drivePath + thisWorkPath + "/";
  } else {
    workPath = homePath + thisWorkPath;
    externalWorkPath = drivePath + thisWorkPath;
  }

  rootPath = workPath + "rootFiles/" + dataSubDir;
  dataPath = externalWorkPath + "data/" + dataSubDir;
  plotPath = workPath + "Plots/" + dataSubDir;
  ptPath = rootPath + "ptData/";
  trigPath = rootPath + "trigData/";
  effPath = rootPath + "effData/";
  xPath = rootPath + "xData/";
  RpPbPath = rootPath + "RpPbData/";

  return;
}


/**
 * Returns a linearly spaced array. The 0th element is lo, and the num-th element is hi.
 */
static double* linspace(double lo, double hi, int num) {
  double* arr = new double[num+1];
  double delta = ((double)(hi)-(double)(lo))/(double)(num);
  for (int i = 0; i <= num; i++) {
    arr[i] = lo + i * delta;
  }
  return arr;
}


/**
 * Returns a logarithmically spaced array, where the 0th element is lo and the num-th element is hi.
 */
static double* logspace(double lo, double hi, int num) {
  double loghi = TMath::Log2(hi);
  if (lo == 0) {
    double* arr = linspace(TMath::Log2(hi/(100*num)), loghi, num);
    for (int i = 0; i <= num; i++) {
      arr[i] = TMath::Power(2, arr[i]);
    }
    return arr;
  } else {
    double loglo = TMath::Log2(lo);
    double* arr = linspace(loglo, loghi, num);
    for (int i = 0; i <= num; i++) {
      arr[i] = TMath::Power(2, arr[i]);
    }
    return arr;
  }
}
