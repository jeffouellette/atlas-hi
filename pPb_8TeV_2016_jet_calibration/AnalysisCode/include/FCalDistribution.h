#ifndef __FCalDistribution_h__
#define __FCalDistribution_h__

#include <TString.h>

namespace pPb8TeV2016JetCalibration {

void FCalDistribution (const int dataSet,
                       const double luminosity = 0, 
                       const bool isMC = false,
                       const bool isPeriodA = false, 
                       const TString inFileName = "",
                       const double crossSection_microbarns = 0,
                       const double filterEfficiency = 0,
                       const int numberEvents = 0);

} // end namespace

#endif
