#ifndef __Intercalibration_h__
#define __Intercalibration_h__

namespace JetCalibration {

/**
 * Primary macro.
 * dataSet: Data set identifier. This should be a run number for data or some other identifier for MC (e.g., slice number).
 * isPeriodA: flag that is raised for MC (meaningless if isMC is false)
 * inFileName: Input root file name where tree is stored; if == "" code will try to guess file name based on other info
 * crossSection_microbarns: Total cross section of the process. Should only be defined for MC.
 * filterEfficiency: Filtering efficiency in the MC generation.
 * numberEvents: The total number of events generated in the MC sample.
 */
void Intercalibration (const char* directory,
                       const int dataSet,
                       const bool isPeriodA,
                       const double pt_low,
                       const double pt_high,
                       const char* inFileName,
                       const bool cutTruthEgammas,
                       const double crossSection_microbarns = 0,
                       const double filterEfficiency = 0,
                       const int numberEvents = 0);

} // end namespace

#endif