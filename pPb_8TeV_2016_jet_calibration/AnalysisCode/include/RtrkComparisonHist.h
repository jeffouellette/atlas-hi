#ifndef __RtrkComparisonHist_h__
#define __RtrkComparisonHist_h__

#include "Params.h"
#include <Initialization.h>
#include <TH2D.h>
#include <TString.h>

using namespace atlashi;

namespace pPb8TeV2016JetCalibration {

TH1D* GetProfileX(const TString name, TH2D* hist, const int nbinsx, const double* xbins, const bool useFit);

TH1D* GetDataOverMC(const TString name, TH2D* data, TH2D* mc, const int numxbins, const double* xbins, const int numybins, const double* ybins, const bool useFit);

void RtrkComparisonHist ();

} // end namespace

#endif
