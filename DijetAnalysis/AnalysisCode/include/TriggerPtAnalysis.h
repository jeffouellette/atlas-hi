#ifndef __TriggerPtAnalysis_h__
#define __TriggerPtAnalysis_h__

namespace JetAnalysis {

// IDEA: plot the number of times each trigger fired for a particular run number.
// For overlapping triggers, this will improve statistics by choosing the trigger that fires more often.
void TriggerPtAnalysis(const int thisRunNumber, // Run number identifier.
                       double luminosity); // Integrated luminosity for this run. Presumed constant over the run period.

} // end namespace

#endif
