#ifndef __Utils_h__
#define __Utils_h__

#include <GlobalParams.h>


namespace pPb8TeV2016JetCalibration {

TH1D* GetProfileX (const TString name, TH2D* hist, const int nbinsx, const double* xbins, const bool useFit = false);

TH1D* GetProfileY (const TString name, TH2D* hist, const int nbinsy, const double* ybins, const bool useFit = false);

//TH1D* GetDataOverMC (const TString name, TH2D* data, TH2D* mc, const int numxbins, const double* xbins, const int numybins, const double* ybins, const bool useFit);
TH1D* GetDataOverMC (const TString name, TH2D* data, TH2D* mc, const int numbins, const double* bins, const bool useFit = false, const TString axis = "x");

} // end namespace

#endif
