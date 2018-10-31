#ifndef __TriggerEfficiencyAnalysisHist_h__
#define __TriggerEfficiencyAnalysisHist_h__

namespace JetAnalysis {

/**
 * Bootstraps the trigger efficiency curve in pt.
 */
void bootstrap(Trigger* trig, TGraphAsymmErrors** effGraphArr);


void TriggerEfficiencyAnalysisHist();

} // end namespace

#endif
