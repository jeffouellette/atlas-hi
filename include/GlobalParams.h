#ifndef __GlobalParams_h__
#define __GlobalParams_h__

#include <TMath.h>
#include <TString.h>

using namespace std;

/**
 * Contains useful variables and directory information for 2016 pPb data analyses.
 * Author: Jeff Ouellette
 * Dated: 4/25/2018
 */

/** User defined parameters **/

const bool debugStatements = false; // Print out periodic statements to monitor code flow
const TString homePath = std::getenv ("ATLAS_PATH"); // ATLAS Heavy Ions home directory (stores code, pdfs, etc.)
const TString drivePath = std::getenv ("ATLAS_DATA_PATH"); // ATLAS Heavy Ions external drive directory (stores grid output)
const TString intPath = std::getenv ("ATLAS_INTERMEDIATE_PATH"); // ATLAS Heavy Ions intermediate files directory (stores smaller root files)

const double dR_HEC = 0.0; // details on the hadronic end cap data cuts.
const double lowerPhiCut = TMath::Pi()-dR_HEC;
const double upperPhiCut = 3.*TMath::Pi()/2.+dR_HEC;
const double lowerEtaCut = 1.5-dR_HEC;
const double upperEtaCut = 3.2+dR_HEC;

const int full_run_list[30] = {313063, 313067, 313100, 313107, 313136, 313187, 313259, 313285, 313295, 313333, 313435, 313572, 313574, 313575, 313603, 313629, 313630, 313688, 313695, 313833, 313878, 313929, 313935, 313984, 314014, 314077, 314105, 314112, 314157, 314170};

/** End user defined parameters **/


/** General (non-user defined) paramters **/

//// More directory information - PLEASE DO NOT CHANGE!!! These values are overwritten when calling triggerUtil::initialize().
//extern TString workPath; // Home analysis directory, should be modified in code outside this path structure
//extern TString externalWorkPath; // External drive storage directory, should be modified in code below
//extern TString intWorkPath; // Base directory for intermediate root files, should be modified in code outside this path structure
//extern TString rootPath; // Where analyzed *.root files are stored. Different analysis modules have different subdirectories here.
//extern TString dataPath; // Where the *.root raw data files (from the CERN grid) are stored.
//extern TString plotPath; // Where plots are stored.
//extern TString ptPath; // Where the pt analysis module output is stored.
//extern TString trigPath;  // Where the trigger fire count module output is stored.
//extern TString effPath; // Where the trigger efficiency module output is stored.
//extern TString xPath; // Where the xa/xp module output is stored.
//extern TString RpPbPath; // Where the R_pPb module output is stored.
//
//extern int numtrigs; // Total number of triggers

// Transverse momentum and pseudorapidity binning
const int MAX_PT = 6000; // Maximum transverse momentum allowed
const double MIN_ETA = -4.9; // Minimum detectable pseudorapidity in hadronic calorimeter
const double MAX_ETA = 4.9; // Maximum detectable pseudorapidity

// Useful constants
const double Z = 82;   // value of Z for Pb
const double A = 208;  // value of A for Pb
const double sqrt_s_nn = 8160; // Collision energy in CoM frame (GeV)
const double electron_mass = 0.000511; // mass of the electron in GeV
const double muon_mass = 0.105658; // mass of the muon in GeV
const double Z_mass = 91.2; // mass of the Z in GeV
const double Z_width = 2.4952; // width of the Z peak in GeV

/** End general parameters **/


#endif
