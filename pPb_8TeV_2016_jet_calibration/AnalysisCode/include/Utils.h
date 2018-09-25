#ifndef __Utils_h__
#define __Utils_h__

#include <GlobalParams.h>


namespace pPb8TeV2016JetCalibration {

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
