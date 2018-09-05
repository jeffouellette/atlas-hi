#ifndef __TriggerEfficiencyAnalysisHist_h__
#define __TriggerEfficiencyAnalysisHist_h__

namespace pPb8TeV2016DijetAnalysis {

/**
 * Bootstraps the trigger efficiency curve in pt.
 */
void bootstrap(Trigger* trig, TGraphAsymmErrors** effGraphArr);


void TriggerEfficiencyAnalysisHist();

} // end namespace

#endif
