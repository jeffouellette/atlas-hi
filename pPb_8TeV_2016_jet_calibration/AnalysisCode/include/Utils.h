#ifndef __Utils_h__
#define __Utils_h__

#include <GlobalParams.h>
#include <TFile.h>

namespace pPb8TeV2016JetCalibration {

/**
 * Returns the appropriate file in the given directory.
 * For MC, inFileName MUST be specified.
 */
TFile* GetFile (const char* directory, const int dataSet, const bool isMC, const char* inFileName = "");


/**
 * Returns an abbreviated, unique identifier for a given dataset.
 */
TString GetIdentifier (const int dataSet, const char* inFileName, const bool isMC, const bool isSignalOnlySample, const bool periodA);


/**
 * Returns the initial systematic error on the 2015 cross-calibration as a function of jet pT and eta.
 * Requires xCalibSystematicsFile to be defined and open, else will return 0.
 */
double GetXCalibSystematicError (const double jpt, const double jeta);


/**
 * Returns the additional systematics associated with applying the 2015 cross-calibration to the 2016 pPb, as a function of jet pT and eta.
 * Requires dataOverMCFile to be defined and open, else will return 0.
 */
double GetNewXCalibSystematicError (const double jeta, const double refpt, const bool periodA);


/**
 * Separates each point on a TGraphAsymmErrors by delta along the x axis, so that the errors don't overlap.
 */
void deltaize (TGraphAsymmErrors* tg, const double delta = 0, const bool logx = false);


/**
 * Makes a TGraphAsymmErrors from the input histogram.
 */
TGraphAsymmErrors* make_graph (TH1* h, const float cutoff = -1);


/**
 * Returns the TProfile of an input histogram along the x axis. Can use either statistical mean or gaussian mean.
 */
TH1D* GetProfileX (const TString name, TH2D* hist, const int nbinsx, const double* xbins, const bool useFit = false);


/**
 * Returns the TProfile of an input histogram along the y axis. Can use either statistical mean or gaussian mean.
 */
TH1D* GetProfileY (const TString name, TH2D* hist, const int nbinsy, const double* ybins, const bool useFit = false);


/**
 * Returns a histogram with the TProfile of data over the TProfile of MC along either the x or y axes. Can use either the statistical or gaussian mean.
 */
TH1D* GetDataOverMC (const TString name, TH2D* data, TH2D* mc, const int numbins, const double* bins, const bool useFit = false, const TString axis = "x");

} // end namespace

#endif
