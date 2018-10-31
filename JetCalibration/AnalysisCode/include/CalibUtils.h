#ifndef __CalibUtils_h__
#define __CalibUtils_h__

#include <GlobalParams.h>
#include <TFile.h>

namespace JetCalibration {

extern TFile* xCalibSystematicsFile;
extern TFile* dataOverMCFile;
extern TFile* purityFile;

/**
 * Returns the initial systematic error on the 2015 cross-calibration as a function of jet pT and eta.
 * Requires xCalibSystematicsFile to be defined and open, else will return 0.
 */
double GetXCalibSystematicError (const double jpt, const double jeta);


/**
 * Returns the photon purities as a function of photon pT and eta.
 * Requires purityFile to be defined and open, otherwise will return 1 (perfect purity).
 */
double GetPurity (const double ppt, const double peta);


/**
 * Returns the additional systematics associated with applying the 2015 cross-calibration to the 2016 pPb, as a function of jet pT and eta.
 * Requires dataOverMCFile to be defined and open, else will return 0.
 */
double GetNewXCalibSystematicError (const double jeta, const double refpt, const bool periodA);

} // end namespace

#endif
