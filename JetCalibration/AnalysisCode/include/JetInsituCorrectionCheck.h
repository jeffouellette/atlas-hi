#ifndef __JetInsituCorrectionCheck_h__
#define __JetInsituCorrectionCheck_h__

namespace JetCalibration {

/**
 * Primary macro.
 * dataSet: Data set identifier. This should be a run number for data or some other identifier for MC (e.g., slice number).
 * luminosity: Integrated luminosity for this run. Presumed constant over the run period. Meaningless for MC.
 * isPeriodA: flag that is raised for period A
 */
void JetInsituCorrectionCheck (const char* directory,
                               const int dataSet,
                               const double luminosity = 0, 
                               const bool isPeriodA = false);

} // end namespace

#endif
