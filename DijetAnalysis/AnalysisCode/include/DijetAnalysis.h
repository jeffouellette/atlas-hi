#ifndef __DijetAnalysis_h__
#define __DijetAnalysis_h__

#include <GlobalParams.h>
#include <Utils.h>

using namespace atlashi;

namespace JetAnalysis {

const int numxbins = 50;
const int numqbins = 8;
const int numq2bins = 120;
const int numq2xbins = 100;
const int nummbins = 50;
const int numfcalbins = 60;

const double* xbins = logspace(1.6e-4, 1.6, numxbins);
const double* qbins = logspace(20, 1200, numqbins);
const double* q2bins = logspace(5e2, 3e6, numq2bins);
const double* q2xbins = logspace(6e-5, 2, numq2xbins);
const double* mbins = logspace(20, 2500, nummbins);
const double* fcalbins = logspace(10, 500, numfcalbins);

void DijetAnalysis(const int dataSet, // Data set identifier. If not MC, this should be a run number. If MC, this should be whatever number follows "tid" in the MC file.
                   const double luminosity, // Integrated luminosity for this run. Presumed constant over the run period.
                   const bool isMC, // MC flag
                   const bool isPeriodA); // period A flag

} // end namespace

#endif
