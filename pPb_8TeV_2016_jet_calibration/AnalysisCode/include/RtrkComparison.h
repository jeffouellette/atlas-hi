#ifndef __RtrkComparison_h__
#define __RtrkComparison_h__

#include <TString.h>

namespace pPb8TeV2016JetCalibration {

/**
 * Calculates the original systematic error on this jet from the cross calib.
 * jpt: pt of the jet
 * jeta: eta of the jet
 */
double GetXCalibSystematicError (const double jpt, const double jeta);


/**
 * Primary macro.
 * dataSet: Data set identifier. This should be a run number for data or some other identifier for MC (e.g., slice number).
 * luminosity: Integrated luminosity for this run. Presumed constant over the run period. Meaningless for MC.
 * isMC: is data/MC flag.
 * isPeriodA: flag that is raised for MC (meaningless if isMC is false)
 * inFileName: Input root file name where tree is stored; if == "" code will try to guess file name based on other info
 */
void RtrkComparison (const int dataSet,
                     const double luminosity = 0, 
                     const bool isMC = false,
                     const bool isPeriodA = false, 
                     const TString inFileName = "",
                     const double crossSection_microbarns = 0,
                     const double filterEfficiency = 0,
                     const int numberEvents = 0);

} // end namespace

#endif
