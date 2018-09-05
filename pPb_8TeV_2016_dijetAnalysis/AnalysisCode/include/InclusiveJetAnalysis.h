#ifndef __InclusiveJetAnalysis_h__
#define __InclusiveJetAnalysis_h__

namespace pPb8TeV2016DijetAnalysis {

void InclusiveJetAnalysis(const int dataSet, // Run number identifier.
                          const double luminosity, // Integrated luminosity for this run. Presumed constant over the run period.
                          const bool isMC,
                          const bool isPeriodA);

} // end namespace

#endif
