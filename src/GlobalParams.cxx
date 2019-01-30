#include "GlobalParams.h"

/**
 * Contains useful variables and directory information for 2016 pPb data analyses.
 * Author: Jeff Ouellette
 * Dated: 4/25/2018
 */

namespace atlashi {

TString workPath = ""; // Home analysis directory, should be modified in code outside this path structure
TString externalWorkPath = ""; // External drive storage directory, should be modified in code below
TString rootPath = ""; // Where analyzed *.root files are stored. Different analysis modules have different subdirectories here.
TString dataPath = ""; // Where the *.root raw data files (from the CERN grid) are stored.
TString plotPath = ""; // Where plots are stored.
TString ptPath = ""; // Where the pt analysis module output is stored.
TString trigPath = "";  // Where the trigger fire count module output is stored.
TString effPath = ""; // Where the trigger efficiency module output is stored.
TString xPath = ""; // Where the xa/xp module output is stored.
TString RpPbPath = ""; // Where the R_pPb module output is stored.

int numtrigs = 0; // number of triggers that have been initialized

} // end namespace
