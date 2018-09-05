#ifndef __GlobalParams_h__
#define __GlobalParams_h__

#include <TMath.h>
#include <TH1.h>
#include <TString.h>
#include <TGraphAsymmErrors.h>
#include <TH2D.h>
#include <TH3D.h>

using namespace std;

namespace atlashi {

/**
 * Contains useful variables and directory information for 2016 pPb data analyses.
 * Author: Jeff Ouellette
 * Dated: 4/25/2018
 */

/** User defined parameters **/

extern const bool debugStatements = false; // Print out periodic statements to monitor code flow
extern const TString homePath = "/Users/jeffouellette/Research/atlas-hi/"; // ATLAS Heavy Ions home directory
extern const TString drivePath = "/Volumes/My Passport/Research/atlas-hi/"; // ATLAS Heavy Ions external drive directory

extern const double dR_HEC = 0.4; // details on the hadronic end cap data cuts.
extern const double lowerPhiCut = TMath::Pi()-dR_HEC;
extern const double upperPhiCut = 3.*TMath::Pi()/2.+dR_HEC;
extern const double lowerEtaCut = 1.5-dR_HEC;
extern const double upperEtaCut = 3.2+dR_HEC;

extern const int full_run_list[30] = {313063, 313067, 313100, 313107, 313136, 313187, 313259, 313285, 313295, 313333, 313435, 313572, 313574, 313575, 313603, 313629, 313630, 313688, 313695, 313833, 313878, 313929, 313935, 313984, 314014, 314077, 314105, 314112, 314157, 314170};

/** End user defined parameters **/


/** General (non-user defined) paramters **/

// More directory information - PLEASE DO NOT CHANGE!!! These values are overwritten when calling triggerUtil::initialize().
extern TString workPath; // Home analysis directory, should be modified in code outside this path structure
extern TString externalWorkPath; // External drive storage directory, should be modified in code below
extern TString rootPath; // Where analyzed *.root files are stored. Different analysis modules have different subdirectories here.
extern TString dataPath; // Where the *.root raw data files (from the CERN grid) are stored.
extern TString plotPath; // Where plots are stored.
extern TString ptPath; // Where the pt analysis module output is stored.
extern TString trigPath;  // Where the trigger fire count module output is stored.
extern TString effPath; // Where the trigger efficiency module output is stored.
extern TString xPath; // Where the xa/xp module output is stored.
extern TString RpPbPath; // Where the R_pPb module output is stored.

// Transverse momentum and pseudorapidity binning
extern const int MAX_PT = 6000; // Maximum transverse momentum
extern const double MIN_ETA = -4.9; // Minimum detectable pseudorapidity in hadronic calorimeter
extern const double MAX_ETA = 4.9; // Maximum detectable pseudorapidity

extern int numtrigs; // Total number of triggers

// Useful constants
extern const double Z = 82;   // value of Z for Pb
extern const double A = 208;  // value of A for Pb
extern const double sqrt_s_nn = 8160; // Collision energy in CoM frame (GeV)
extern const double electron_mass = 0.000511; // mass of the electron in GeV
extern const double muon_mass = 0.105658; // mass of the muon in GeV
extern const double Z_mass = 91.2; // mass of the Z in GeV
extern const double Z_width = 2.4952; // width of the Z peak in GeV

extern const double pi = TMath::Pi();

/** End general parameters **/


/**
 * Returns a TString summarizing a measurement.
 * By default, there is only 1 significant digit in the error.
 * E.g., FormatMeasurement (40.58, 1.29) returns "40#pm1".
 * Or, FormatMeasurement (40.58, 1.29, 2) returns "40.6#pm1.3".
 */
const char* FormatMeasurement (double val, double err, const int n=1);


/**
 * Modifies the directory strings to point to the correct locations.
 */
void SetupDirectories (const TString dataSubDir, const TString thisWorkPath);


/**
 * Returns a linearly spaced array. The 0th element is lo, and the num-th element is hi.
 */
double* linspace (double lo, double hi, int num);


/**
 * Returns a logarithmically spaced array, where the 0th element is lo and the num-th element is hi.
 */
double* logspace (double lo, double hi, int num);


/**
 * Returns the equivalent angle in the range 0 to 2pi.
 */
double InTwoPi (double phi);


/**
 * Returns the difference between two angles in 0 to pi.
 */
double DeltaPhi (double phi1, double phi2);


/**
 * Returns dR between two eta, phi coordinates.
 */
double DeltaR (const double eta1, const double eta2, const double phi1, const double phi2 );


/**
 * Returns true iff this eta, phi coordinate lies in the disabled HEC region.
 */
bool InDisabledHEC (const double eta, double phi);


/**
 * Returns true iff this eta lies within the EMCal.
 */
bool InEMCal (const float eta);


/**
 * Returns true iff this object is within a given radius in the HCal.
 */
bool InHadCal (const float eta, const float R = 0.4);


/**
 * Calculates the systematic errors on optimal, storing the results in graph.
 */
void CalcSystematics (TGraphAsymmErrors* graph, TH1* optimal, TH1* sys_hi, TH1* sys_lo);

/**
 * Creates a projection of a TH3 with the specified axes by integrating between min and max on the 3rd axis of the TH3.
 */
TH2D* Project2D (TString name, TH3D* h3, const TString xaxis, const TString yaxis, const int min, const int max);


} // end namespace

#endif
